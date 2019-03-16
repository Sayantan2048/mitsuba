cd dist
mitsuba -Des=0 -Dbs=16000 -Dabs=0 D:\projects\mitsuba\scenes\veach\veachDirectCv.xml -p 2
mtsgui D:\projects\mitsuba\scenes\veach\veachDirectCv.exr
cd ..