
## This scripts renderes all obj files to respective exr images and then to png from which 
## an mp4 movie is made

# Call mitsuba for each obj file from the folder "objfiles" and save the output into the
# folder "exrfiles"
for i in `seq 0 175`;
do
  mitsuba -Dfilename=./objfiles/levelset_$i.obj -o ./exrfiles/levelset_$(printf %03d $i).exr test.xml
done

# Using 'exrtopng' application in Ubuntu, all images are converted into png and stored in
# the folder "pngfiles"
for i in `seq 0 175`;
do
  exrtopng ./exrfiles/levelset_$(printf %03d $i).exr ./pngfiles/levelset_$(printf %03d $i).png
done

# Using ffmpeg make a movie
cd pngfiles
ffmpeg -qscale 5 -r 20 -b 9600 -i levelset_%03d.png crystal.mp4

