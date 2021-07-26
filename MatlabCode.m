%Embedding process to get WaterMarked Image

%Read Hospital Logo image. Convert it to grayscale. Get binary Vector W1%
logo = imread('logo.jpeg');
%[file, path] = uigetfile({'.'},'File Selector');
%img_path = strcat(path, file);
%logo = im2double(imread(img_path));
logo = rgb2gray(logo);
logo = imresize(logo,[64 64]);
%disp(size(logo));

b_img = imbinarize(logo);
W1 = b_img(:);


%Read Patient Info file. Convert characters to ascii values and then get binary array W2 of ascii values.%

file = fileread('p_info.txt');


ascii_vals = double(file);


ascii_array = cell2mat(arrayfun(@(c) bitget(c, 8:-1:1), ascii_vals, 'UniformOutput', false));
W2 = ascii_array(:);

%Setting LSB's to 0 and computing hash funciton of image.%
%[file, path] = uigetfile({'.'},'File Selector');
%img_path = strcat(path, file);
img = im2double(imread('project.jpeg'));
img = rgb2gray(img);

lsb_img = dec2bin(img);
lsb_img(:,8) = '0';
lsb_img = bin2dec(lsb_img);
lsb_img = reshape(lsb_img, size(img));

hash = DataHash(lsb_img);


%Convert hash to binary vector%
hash_ascii = double(hash);
hash_ascii_array = cell2mat(arrayfun(@(c) bitget(c, 8:-1:1), hash_ascii, 'UniformOutput', false));
W3 = hash_ascii_array(:);


%Concatenate all three%
W = cat(1, W1, W2, W3);



%W* = W xor P%

rng(10);
d = size(W);
P = rand(d(1),d(2)) < 0.5;
W_star = xor(W,P);



%%%%%%%%%%%%%%%%%%%%%% ROI and RONI


[file, path] = uigetfile({'.'},'File Selector');
img_path = strcat(path, file);
original = im2double(imread(img_path));
I = im2double(imread(img_path));
buffer = im2double(imread(img_path));
I = rgb2gray(I); %converting 3D to 2D 
%figure, imshow(I);
Tfinal = graythresh(buffer); %Where Tfinal is the threshold Value

[rows,cols] = size(original);
%Finding binarized image by Otsue Method Calculated Tfinal Value.
%Step (4)
border = 4; % I am taking a border of 4
for i=border+1:rows-border
    for j=border+1:cols-border
        if(original(i,j) > Tfinal)
            original(i,j) = 255; 
        end
    end
end 
original = rgb2gray(original);

[m,n] = size(original);
disp(m);
disp(n);
x=200; y=259;
Right_lung = regiongrowing(original,x,y,0.2);
x=200; y=100;
Left_lung = regiongrowing(original,x,y,0.2);



tagged_img = zeros(size(original));
tagged_img = imbinarize(tagged_img,Tfinal);
for row = 1:m
    for col = 1:n
        if (Right_lung(row,col) == 1)
            tagged_img(row,col) = 1;
        end
    end
end
for row = 1:m
    for col = 1:n
        if (Left_lung(row,col) == 1)
            tagged_img(row,col) = 1;
        end
    end
end

%figure,imshow(tagged_img);
ROI = original;

for row = 1:m
    for col = 1:n
        if (tagged_img(row,col) == 0)
            ROI(row,col) = 255;
        end
    end
end
original = im2double(imread(img_path));
original = rgb2gray(original);
RONI = original;
title('ROI');
figure,imshow(ROI);



for row = 1:m
    for col = 1:n
        if (ROI(row,col) ~= 255)
            RONI(row,col) = 255;
        end
    end
end
title('ROI');
figure,imshow(RONI);

disp(size(RONI));

rng(10);
im = RONI;

vec = randperm(numel(im));
vec = reshape(vec, size(im));
out = im(vec);

[W_r,W_c] = size(W_star);
%out_bin = dec2bin(out);
 
disp('out');
disp(size(out));
disp(W_r);
% disp(size(W_star));
% disp(size(out_bin));
% disp('M');
% disp(m);

itr_col = 1;
itr_row = 1;
for row = 1:W_r
     if (itr_row == m)
        itr_row = 1;
        itr_col = itr_col+1;
        disp('yes');
     end
     temp = out(itr_row,itr_col);
%      disp('temp');
%      disp(size(temp));
%      %temp_bin = de2bi(temp);
%      temp_bin = myDec2Bin(temp);
%      disp('temp_bin');
%      disp(temp_bin);
%      temp_bin(1,8) =  W_star(row,1);
%      temp_dec = myBin2Dec(temp_bin);
%      temp_dec = reshape(temp_dec,size(temp));
     temp_bin = dec2bin(temp);
     [trow,tcolumn]=size(temp_bin);
     %temp_bin(tcolumn)=num2str(W_star(row,1));
     if(W_star(row,1) == 0)
         temp_bin(tcolumn)= '0';
     else
         temp_bin(tcolumn)= '1';
     end
     disp(temp_bin);
     temp_dec = bin2dec(temp_bin);
     
     out(itr_row,itr_col) = temp_dec;
     itr_row = itr_row + 1;
     %disp(size(temp_bin));
     
 end
% 

%To reconstruct the original image
reconstruct = zeros(size(im), class(im));
reconstruct(vec) = out;

title('RONI');
figure,imshow(out);
title('scrambled RONI');
figure,imshow(reconstruct);

%reuniting ROI and RONI

for row = 1:m
    for col = 1:n
        if (ROI(row,col) ~= 255)
            reconstruct(row,col) = ROI(row,col);
        end
    end
end
title('unscrambled RONI');
figure,imshow(reconstruct);

imwrite(reconstruct, 'D:/temp1.tif');
WaterMarkedImage = imread('D:/temp1.tif');
title('Water Marked Image');
figure,imshow(WaterMarkedImage);



%%%%Extraction Process



[file, path] = uigetfile({'.'},'File Selector');
img_path = strcat(path, file);
original = im2double(imread(img_path));
I = im2double(imread(img_path));
buffer = im2double(imread(img_path));

%figure, imshow(I);
Tfinal = graythresh(buffer); %Where Tfinal is the threshold Value

[rows,cols] = size(original);
%Finding binarized image by Otsue Method Calculated Tfinal Value.
%Step (4)
border = 4; % I am taking a border of 4
for i=border+1:rows-border
    for j=border+1:cols-border
        if(original(i,j) > Tfinal)
            original(i,j) = 255; 
        end
    end
end 

[m,n] = size(original);
disp(m);
disp(n);
x=200; y=259;
Right_lung = regiongrowing(original,x,y,0.2);
x=200; y=100;
Left_lung = regiongrowing(original,x,y,0.2);



tagged_img = zeros(size(original));
tagged_img = imbinarize(tagged_img,Tfinal);
for row = 1:m
    for col = 1:n
        if (Right_lung(row,col) == 1)
            tagged_img(row,col) = 1;
        end
    end
end
for row = 1:m
    for col = 1:n
        if (Left_lung(row,col) == 1)
            tagged_img(row,col) = 1;
        end
    end
end

%figure,imshow(tagged_img);
ROI = original;

for row = 1:m
    for col = 1:n
        if (tagged_img(row,col) == 0)
            ROI(row,col) = 255;
        end
    end
end
original = im2double(imread(img_path));
RONI = original;


for row = 1:m
    for col = 1:n
        if (ROI(row,col) ~= 255)
            RONI(row,col) = 255;
        end
    end
end


rng(10);
im = RONI;

vec = randperm(numel(im));
vec = reshape(vec, size(im));
out = im(vec);

W_star_extraction = zeros(size(W_star));
[W_r,W_c] = size(W_star);
itr_col = 1;
itr_row = 1;

for row = 1:W_r
     if (itr_row == m)
        itr_row = 1;
        itr_col = itr_col+1;
        disp('yes');
     end
     temp = out(itr_row,itr_col);

     temp_bin = dec2bin(temp);
     [trow,tcolumn]=size(temp_bin);
     %temp_bin(tcolumn)=num2str(W_star(row,1));
     if(temp_bin(tcolumn) == '0')
          W_star(row,1) = 0;
     else
         W_star(row,1) = 1;
     end
     
     temp_dec = bin2dec(temp_bin);
     
     out(itr_row,itr_col) = temp_dec;
     itr_row = itr_row + 1;
end
%The extracted LSBS from the Water Marked Image.
disp(W_star);
%To reconstruct the original image
reconstruct = zeros(size(im), class(im));
reconstruct(vec) = out;



