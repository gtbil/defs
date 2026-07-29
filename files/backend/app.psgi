# Plack replacement for the old Apache2 + mod_perl setup.
#
#   /cgi-bin/*   -> executed as CGI (Plack::App::CGIBin)
#   /            -> static files from the htdocs tree (results + assets)
#
# Directories come from BP3_WEB (the writable copy staged by the Electron
# shell). This removes the need for Apache, mod_perl, sites-available,
# ports.conf, and running anything as root.

use strict;
use warnings;
use Plack::Builder;
use Plack::App::CGIBin;
use Plack::App::File;

my $web = $ENV{BP3_WEB} or die "BP3_WEB not set";
my $cgi_dir    = "$web/cgi-bin";
my $htdocs_dir = "$web/htdocs";

# These legacy scripts rely on old-style "one process per request" CGI
# semantics (top-level lexicals closed over by named subs like
# input_screen()). Plack::App::CGIBin's default in-process mode compiles
# each script once and reuses it across requests, which breaks that
# closure pattern after the first request (primer defaults/URLs render
# blank). Force real fork+exec per request instead.
my $cgi_app = Plack::App::CGIBin->new(
    root    => $cgi_dir,
    exec_cb => sub { 1 },
)->to_app;
my $static  = Plack::App::File->new(root => $htdocs_dir)->to_app;

builder {
    mount "/cgi-bin" => $cgi_app;
    mount "/"        => $static;
};
