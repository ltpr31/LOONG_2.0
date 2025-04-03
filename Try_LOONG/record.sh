
###install LOONG
sudo dpkg -i /path/to/LOONG_V2.0.deb
echo 'export LD_LIBRARY_PATH=/usr/local/lib/LOONG/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
echo 'export PATH=$PATH:/usr/local/LOONG/bin' >> ~/.bashrc
source ~/.bashrc



###usage command
draw_a_loong_protect  --rootPath  /path/to/try_LOONG_here  \
--txtName A1_S1_tissue_positions_list.txt  \
--imgName tissue_hires_image.png  \
--outTxtPrefix albst \
--row_std 10 --outImg_height 500  \
--x_scale_factor 0.2019998  --y_scale_factor 0.2019998  \
--point_size 4  --line_size 2  --splitImg splitted.png


##Multiple layers has been divided through "--splitImg splitted.png".
##If you want to divide them again, just delete "--splitImg splitted.png"



