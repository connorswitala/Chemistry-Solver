# Build nomenclature from .nlo → .nls
add_cus_dep('nlo','nls',0,'nomencl');
sub nomencl {
  my ($base) = @_;
  system("makeindex -s nomencl.ist -o '$base.nls' '$base.nlo'");
}