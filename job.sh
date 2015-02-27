#!/bin/sh
D=/afs/cern.ch/user/i/inaki
cp $D/LiE/Lie.exe .
cp $D/LiE/lie .
cp $D/LiE/symtensor21.in .
./lie < symtensor21.in > symtensor21.out.batch
cp symtensor21.out.batch $D/batch_result