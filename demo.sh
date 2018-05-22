#!/bin/bash
./memoryfill 15G
./copmem -l 100 -o hm-100.txt hum.all.fa mus.all.fa
./memoryfill 15G
./copmem -l 200 -o hm-200.txt hum.all.fa mus.all.fa
./memoryfill 15G
./copmem -l 300 -o hm-300.txt hum.all.fa mus.all.fa
./memoryfill 15G
./copmem -l 400 -o hm-400.txt hum.all.fa mus.all.fa
