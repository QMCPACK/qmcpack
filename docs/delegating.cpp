/** Generic delgate input class.
 *  The delegate input requires:
 *  A constructor that takes an xmlNodePtr
 *  A factoryFunction that matches the signature discussed in InputSection::DelegateHandler
 */

class AnotherInput
{
public:
  class AnotherInputSection : public InputSection
  {
  ...
  };

  AnotherInput(xmlNodePtr cur)
  {
  ...
  }
};

/// Factory function for the delegate Input class.
std::any makeAnotherInput(xmlNodePtr cur, std::string& value_label)
{
  AnotherInput another_input{cur};
  value_label = another_input.get_name();
  return another_input;
}

/** Input class delegating to other input classes.
 */
class DelegatingInput
{
  class DelegatingInputSection : public InputSection
  {
  public:
    DelegatingInputSection()
    {
      section_name = "DelegateTest";
      ...
      delegates    = {"anotherinput"};
      InputSection::registerDelegate("anotherinput", makeAnotherInput);
    }
  };

public:
  DelegatingInput(xmlNodePtr cur) { dins_.readXML(cur); }

  // rather than write a real input class we just pass through InputSection for this test.
  bool has(const std::string& name) const { return dins_.has(name); }
  template<typename T>
  T get(const std::string& name) const
  {
    return dins_.get<T>(name);
  }

private:
  std::string name_ = "DelegatingInput";
  DelegatingInputSection dins_;
};
