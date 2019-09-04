#!/bin/bash

# Run all the cases in one go

echo "----------"
echo "f01oneline"
echo "----------"
../miniRT.py --atm ../inputs/oneline.atm --par ../inputs/oneline.par --fname oneline --wnlo 4264.0 --wnhi 4455.0 --wnsamp 0.5 --osamp 2160 --outdir ../output/f01oneline/

echo "----------"
echo "f02fewline"
echo "----------"
../miniRT.py --atm ../inputs/fewline.atm --par ../inputs/fewline.par --fname fewline --wnlo 3030.0 --wnhi 5000.0 --wnsamp 0.5 --osamp 2160 --outdir ../output/f02fewline/

echo "------------"
echo "f03multiline"
echo "------------"
../miniRT.py --atm ../inputs/multiline.atm --par ../inputs/multiline.par --fname multiline --wnlo 2500.0 --wnhi 5000.0 --wnsamp 0.25 --osamp 2160 --outdir ../output/f03multiline/

echo "-------------"
echo "f04broadening"
echo "-------------"
../miniRT.py --atm ../inputs/broadening.atm --par ../inputs/oneline.par --fname broadening --wnlo 4366 --wnhi 4370 --wnsamp 0.005 --osamp 1080 --outdir ../output/f04broadening/

echo "------------"
echo "f05abundance"
echo "------------"
../miniRT.py --atm ../inputs/abundance_atm/0_uniform.atm --par ../inputs/abundance.par --fname 0 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

../miniRT.py --atm ../inputs/abundance_atm/1e-4_uniform.atm --par ../inputs/abundance.par --fname 1 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

../miniRT.py --atm ../inputs/abundance_atm/2e-4_uniform.atm --par ../inputs/abundance.par --fname 2 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

../miniRT.py --atm ../inputs/abundance_atm/3e-4_uniform.atm --par ../inputs/abundance.par --fname 3 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

../miniRT.py --atm ../inputs/abundance_atm/4e-4_uniform.atm --par ../inputs/abundance.par --fname 4 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

../miniRT.py --atm ../inputs/abundance_atm/5e-4_uniform.atm --par ../inputs/abundance.par --fname 5 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

../miniRT.py --atm ../inputs/abundance_atm/6e-4_uniform.atm --par ../inputs/abundance.par --fname 6 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

../miniRT.py --atm ../inputs/abundance_atm/7e-4_uniform.atm --par ../inputs/abundance.par --fname 7 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

../miniRT.py --atm ../inputs/abundance_atm/8e-4_uniform.atm --par ../inputs/abundance.par --fname 8 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

../miniRT.py --atm ../inputs/abundance_atm/9e-4_uniform.atm --par ../inputs/abundance.par --fname 9 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

../miniRT.py --atm ../inputs/abundance_atm/1e-3_uniform.atm --par ../inputs/abundance.par --fname 10 --wnlo 4365 --wnhi 4371 --wnsamp 0.005 --osamp 1080 --outdir ../output/f05abundance/

echo "-----------"
echo "f06blending"
echo "-----------"
../miniRT.py --atm ../inputs/blending.atm --par ../inputs/blending.par --fname blending --wnlo 4366 --wnhi 4370 --wnsamp 0.005 --osamp 1080 --outdir ../output/f06blending/
