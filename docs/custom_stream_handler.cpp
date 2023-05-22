class CustomTypeInput : public InputSection
{
public:
  struct LabeledPosition
  {
    std::string letters;
    std::array<int, 3> numbers;
  };
  CustomTestInput()
  {
    ...
    attributes   = {"name", "custom_attribute"};
    custom       = {"labeled_position", "custom_attribute"};
  }
  void setFromStreamCustom(const std::string& ename, const std::string& name, std::istringstream& svalue) override
  {
    if (ename == "labeled_position")
    {
      LabeledPosition ws;
      svalue >> ws.letters;
      svalue >> ws.numbers[0];
      svalue >> ws.numbers[1];
      svalue >> ws.numbers[2];
      values_[name] = ws;
    }
    else if (name == "custom_attribute")
    {
      std::string cus_at;
      // otherwise you will get not consume the entire attribute just the next piece as determine by operator>> and
      // the type you are going into.
      std::getline(svalue, cus_at);
      ... // do your custom parsing.
      values_[name] = cus_at;
    }
    else
      throw std::runtime_error("bad name passed: " + name +
                               " or custom setFromStream not implemented in derived class.");
  }
};
