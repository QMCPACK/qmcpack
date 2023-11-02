#!/bin/bash

PPCONVERT=$1

${PPCONVERT} --gamess_pot O.BFD.gamess --s_ref "1s(2)2p(4)"  --p_ref "1s(2)2p(4)" --d_ref "1s(2)2p(4)" --xml O.BFD.xml

exit $?
# if (( $? != 0 )) then
#    echo "Fail"
# fi
