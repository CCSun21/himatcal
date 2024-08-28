#!/bin/bash
python -m pyqrc EC_TS-c1s2.log --nproc 32 --mem 64GB --freqnum 1 --amp -0.3 --route 'B3LYP/6-31G*' --name QRCR
python -m pyqrc EC_TS-c1s2.log --nproc 32 --mem 64GB --freqnum 1 --amp 0.3 --route 'B3LYP/6-31G*' --name QRCF
