import testing
import os
from wannier90_input import Wannier90Input

associated_files = dict()

input_files = [
    'BaTiO3.win',
    'Fe.win',
    'LaVO3.win',
    'Pt.win',
    'Si.win',
    'W.win',
    'benzene.win',
    'cnt55.win',
    'copper.win',
    'diamond.win',
    'gaas.win',
    'graphite.win',
    'iron_dn.win',
    'iron_up.win',
    'lead.win',
    'silane.win',
    'silicon.win'
]

def get_files():
    return testing.collect_unit_test_file_paths('wannier90_input',associated_files)
#end def get_files


def test_files():
    filenames = input_files
    files = get_files()
    assert(set(files.keys())==set(filenames))
#end def test_files

def pattern_in_text(pattern, text):
    """Check if the pattern exists in the text."""
    if pattern.startswith('#') or pattern.startswith('!'):
        # Ignore comment lines
        return True
    else:
        pl = ' '.join(pattern.lower().strip().split()) 
        tl = ' '.join(text.lower().strip().split())
        if pl in tl:
            return True
        elif "=" in pl:
            pattern_key = pl.split('=')[0].strip()
            pattern_value = ' '.join(pl.split('=')[1].strip().split())
            pl_alt = pattern_key + ' : ' + pattern_value
            pl = pattern_key + ' = ' + pattern_value
            return pl in tl or pl_alt in tl
        elif ":" in pl:
            pattern_key = pl.split(':')[0].strip()
            pattern_value = ' '.join(pl.split(':')[1].strip().split())
            pl_alt = pattern_key + ' = ' + pattern_value
            pl = pattern_key + ' : ' + pattern_value
            return pl in tl or pl_alt in tl
        else:
            return False
        
def test_wannier90_input():
    # Directory containing the .win files
    win_files = get_files()
    
    for win_file, win_path in win_files.items():
        # Read the .win file
        with open(win_path, 'r') as f:
            original_content = f.read()
        
        # Create a Wannier90Input object and read the file
        win_input = Wannier90Input()
        win_input.read(win_path)
        
        # Write the input back to a string
        generated_content = ' '.join(win_input.write_text().split())
        expected_patterns = [x.strip() for x in original_content.split('\n') if x.strip()]
        # Compare the generated content to the original
        for pattern in expected_patterns:
            try:
                assert pattern_in_text(pattern, generated_content), f"Missing or incorrect pattern: {pattern}"
            except AssertionError as e:
                print(f"Error for {win_file}: {e}")
                raise e
         
        # print(f"Test passed for {win_file}")

# Run the test
# test_wannier90_input()
