#/bin/bash

cnt=400;
head="<?xml version=\"1.0\"?><VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\"> <Collection>";
tail="</Collection></VTKFile>";
time=0;
step=0.050;

printf "%s\n" $head;

for ((i=0; i<$cnt; i++)) do
    printf "<DataSet timestep=\"%g\" group=\"\" part=\"0\" file=\"/export/home/urawat/Desktop/FTLE/Time_Dependent/DoubleGyre/DoubleGyreA_%.1d.vti\"/>\n" $time $i;
    time=`echo $time + $step | bc`
done

    
printf "%s\n" $tail;
