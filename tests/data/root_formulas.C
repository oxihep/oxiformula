// ROOT's answers for the formulas of root_formulas.txt, for the test in
// ../root_golden.rs. Input: root_formulas.in.tsv (`formula <TAB> to_root`, as
// oxifit-formula writes it); output: root_formulas.tsv. To regenerate, from
// the repository root:
//
//   cargo test -p oxifit-formula --test root_golden -- --ignored write_inputs
//   cd crates/oxifit-formula/tests/data && root -l -b -q root_formulas.C
//
// Each formula is evaluated at x, y, z, t = 1.3, 0.7, 0.2, 0.9 with parameters
// 1, 1.5, 2, 2.5, … by index; then the formula oxifit-formula writes for it
// (to_root) is compiled by ROOT and evaluated the same way.
void root_formulas() {
   gErrorIgnoreLevel = kFatal;
   std::ifstream in("root_formulas.in.tsv");
   std::ofstream out("root_formulas.tsv");
   out << "# formula\tndim\tnpar\tnames\tvalue\tto_root\tto_root_value  (ROOT " << gROOT->GetVersion() << ")\n";
   std::string line;
   double x[4] = {1.3, 0.7, 0.2, 0.9};
   auto eval = [&](TFormula &f) {
      for (int i = 0; i < f.GetNpar(); i++) f.SetParameter(i, 1.0 + 0.5 * i);
      return f.EvalPar(x, nullptr);
   };
   while (std::getline(in, line)) {
      auto tab = line.find('\t');
      std::string src = line.substr(0, tab), root = line.substr(tab + 1);
      TFormula f("f", src.c_str());
      if (!f.IsValid()) { out << src << "\tINVALID\n"; continue; }
      std::string names;
      for (int i = 0; i < f.GetNpar(); i++) { if (i) names += ","; names += f.GetParName(i); }
      double v = eval(f);
      TFormula g("g", root.c_str());
      out << src << "\t" << f.GetNdim() << "\t" << f.GetNpar() << "\t" << names << "\t"
          << TString::Format("%.17g", v) << "\t" << root << "\t"
          << (g.IsValid() ? TString::Format("%.17g", eval(g)) : TString("INVALID")) << "\n";
   }
}
