#!/bin/bash

ppconvert --gamess_pot H.BFD.gamess --s_ref "1s(1)" --log_grid --xml H.xml
ppconvert --gamess_pot O.BFD.CCT.gamess --s_ref "1s(2)2p(4)" --p_ref "1s(2)2p(4)" --d_ref "1s(2)2p(4)" --log_grid --xml O.xml
