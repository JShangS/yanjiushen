clc
clear 
close all

Pc=importdata('d:\Pc.txt');
figure
imagesc(Pc)

Gv_gpu=importdata('d:\GV.txt');
figure
mesh(Gv_gpu)

Gv_cpu=importdata('d:\GV_cpu.txt');
figure
mesh(Gv_cpu)

figure()
mesh(abs(Gv_cpu-Gv_gpu))