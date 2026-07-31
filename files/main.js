// BatchPrimer3 Electron shell
// ---------------------------------------------------------------------------
// Responsibilities:
//   1. Locate the bundled backend (packed conda env + patched CGI web tree).
//   2. On first launch, copy the read-only web tree into a WRITABLE user dir
//      (the legacy CGI writes result files next to itself) and run
//      `conda-unpack` once so the relocated env fixes its own prefixes.
//   3. Pick a localhost port, spawn plackup, poll until it answers.
//   4. Open a BrowserWindow pointed at the CGI entry point.
//   5. Tear the backend child process down on quit.
//
// The same Linux build runs natively on Linux and inside WSL (WSLg) on Windows.

const { app, BrowserWindow, dialog } = require("electron");
const { spawn } = require("child_process");
const http = require("http");
const net = require("net");
const path = require("path");
const fs = require("fs");
const fsp = require("fs/promises");

const PREFERRED_PORT = 8181; // legacy scripts historically assumed this
const START_TIMEOUT_MS = 60_000;

let backendProc = null;
let mainWindow = null;

// --- Resource paths: differ between `npm start` (dev) and a packaged app ----
function resourcesRoot() {
  // In a packaged build, extraResources land in process.resourcesPath.
  // In dev, everything sits next to this file.
  return app.isPackaged ? process.resourcesPath : __dirname;
}

const BUNDLE_DIR = () => path.join(resourcesRoot(), "backend", "bundle"); // packed conda env
const PSGI_FILE = () => path.join(resourcesRoot(), "backend", "app.psgi");
const START_SH = () => path.join(resourcesRoot(), "backend", "start-backend.sh");
const WEB_SRC = () => path.join(resourcesRoot(), "web"); // pristine patched cgi-bin + htdocs

// Writable working copy (legacy CGI needs to write results next to itself)
const WEB_WORK = () => path.join(app.getPath("userData"), "web");

// --- Helpers ----------------------------------------------------------------
function getFreePort(preferred) {
  return new Promise((resolve) => {
    const tryPort = (p, fallback) => {
      const srv = net.createServer();
      srv.once("error", () => (fallback ? tryPort(0, false) : resolve(0)));
      srv.once("listening", () => {
        const { port } = srv.address();
        srv.close(() => resolve(port));
      });
      srv.listen(p, "127.0.0.1");
    };
    tryPort(preferred, true);
  });
}

async function copyDir(src, dst) {
  await fsp.cp(src, dst, { recursive: true });
}

// First-run: stage a writable web tree the CGI can write into.
async function ensureWritableWeb() {
  const work = WEB_WORK();
  const marker = path.join(work, ".staged");
  if (fs.existsSync(marker)) return work;
  await fsp.rm(work, { recursive: true, force: true });
  await copyDir(WEB_SRC(), work);
  await fsp.writeFile(marker, new Date().toISOString());
  return work;
}

function waitForPort(port, deadline) {
  return new Promise((resolve, reject) => {
    const tick = () => {
      const req = http.get(
        { host: "127.0.0.1", port, path: "/", timeout: 2000 },
        (res) => {
          res.destroy();
          resolve();
        }
      );
      req.on("error", () => {
        if (Date.now() > deadline) return reject(new Error("backend did not start in time"));
        setTimeout(tick, 500);
      });
      req.on("timeout", () => {
        req.destroy();
        if (Date.now() > deadline) return reject(new Error("backend timed out"));
        setTimeout(tick, 500);
      });
    };
    tick();
  });
}

async function startBackend() {
  const port = await getFreePort(PREFERRED_PORT);
  const webWork = await ensureWritableWeb();

  const env = {
    ...process.env,
    BP3_PORT: String(port),
    BP3_BUNDLE: BUNDLE_DIR(),
    BP3_PSGI: PSGI_FILE(),
    BP3_WEB: webWork, // writable cgi-bin + htdocs live here
  };

  backendProc = spawn("bash", [START_SH()], {
    env,
    stdio: ["ignore", "pipe", "pipe"],
    detached: true, // own process group, so we can kill the whole tree
  });

  backendProc.stdout.on("data", (d) => process.stdout.write(`[backend] ${d}`));
  backendProc.stderr.on("data", (d) => process.stderr.write(`[backend] ${d}`));
  backendProc.on("exit", (code) => {
    if (code && code !== 0 && mainWindow) {
      dialog.showErrorBox("Backend stopped", `The BatchPrimer3 server exited with code ${code}.`);
    }
  });

  await waitForPort(port, Date.now() + START_TIMEOUT_MS);
  return port;
}

function stopBackend() {
  if (backendProc && !backendProc.killed) {
    try {
      // Negative PID kills the whole process group (plackup + any workers).
      process.kill(-backendProc.pid, "SIGTERM");
    } catch (_) {
      try { backendProc.kill("SIGTERM"); } catch (_) {}
    }
    backendProc = null;
  }
}

// Injected into every loaded CGI page: a floating back button, since the
// plain BrowserWindow has no chrome of its own and the legacy CGI has no
// in-page way to return to a previous step.
const BACK_BUTTON_JS = `
(function () {
  if (document.getElementById("__bp3_back_btn")) return;
  var btn = document.createElement("button");
  btn.id = "__bp3_back_btn";
  btn.textContent = "\\u2039 Back";
  btn.style.cssText = "position:fixed;top:8px;left:8px;z-index:2147483647;" +
    "padding:4px 10px;font:13px sans-serif;cursor:pointer;" +
    "background:#f0f0f0;border:1px solid #999;border-radius:4px;opacity:0.85;";
  btn.disabled = history.length <= 1;
  btn.addEventListener("click", function () { history.back(); });
  document.body.appendChild(btn);
})();
`;

function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1200,
    height: 850,
    title: "BatchPrimer3",
    webPreferences: {
      preload: path.join(__dirname, "preload.js"),
      contextIsolation: true,
      nodeIntegration: false,
    },
  });
  mainWindow.loadFile("loading.html");
  mainWindow.webContents.on("did-finish-load", () => {
    if (mainWindow.webContents.getURL().startsWith("file://")) return; // skip the loading.html splash
    mainWindow.webContents.executeJavaScript(BACK_BUTTON_JS).catch(() => {});
  });
  mainWindow.on("closed", () => (mainWindow = null));
}

app.whenReady().then(async () => {
  createWindow();
  try {
    const port = await startBackend();
    // Entry point of the legacy app. Adjust the path if your patched
    // app.psgi mounts cgi-bin somewhere else.
    await mainWindow.loadURL(`http://127.0.0.1:${port}/cgi-bin/batchprimer3.cgi`);
  } catch (err) {
    dialog.showErrorBox("Could not start BatchPrimer3", String(err && err.message ? err.message : err));
  }
});

app.on("before-quit", stopBackend);
app.on("window-all-closed", () => {
  stopBackend();
  app.quit();
});
process.on("exit", stopBackend);
