#!/env/bin/python
def test_count_kshells():
  from test_physical_system import example_structure_h4
  s1 = example_structure_h4()
  kf = 1.465
  kcut = 5*kf
  nksh = s1.count_kshells(kcut)
  assert nksh == 13
#end def test_count_kshells

if __name__ == '__main__':
  test_count_kshells()
# end __main__
