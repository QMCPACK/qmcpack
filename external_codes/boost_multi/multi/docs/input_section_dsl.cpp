class TestInputSection : public InputSection
{
public:
  TestInputSection()
  {
    section_name   = "Test";
    attributes     = {"name", "samples", "kmax", "full"};
    parameters     = {"label",     "count",   "width",  "rational", "testenum1",
                      "testenum2", "sposets", "center", "density",  "target"};
    required       = {"count", "full"};
    strings        = {"name", "label"};
    multi_strings  = {"sposets"};
    multi_reals    = {"density"};
    multiple       = {"target"};
    integers       = {"samples", "count"};
    reals          = {"kmax", "width"};
    positions      = {"center", "target"};
    bools          = {"full", "rational"};
    enums          = {"testenum1", "testenum2"};
    default_values = {{"name", std::string("demo")},
                      {"samples", int(20)},
                      {"width", Real(1.0)},
                      {"rational", bool(false)}};
  }
}
