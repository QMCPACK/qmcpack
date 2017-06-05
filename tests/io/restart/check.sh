#!/bin/bash
sed "s/#//" qmc_short.s001.scalar.dat | awk '{print $5,$6,$7,$8}' > s001
sed "s/#//" qmc_short.s002.scalar.dat | awk '{print $5,$6,$7,$8}' > s002
diff s001 s002
