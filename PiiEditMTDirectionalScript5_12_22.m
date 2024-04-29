function PiiEditMTDirectionalScript5_12_22 (varargin)
%%
%dispOM3D
%Dmitry Yampolsky 2018

%This code was edited by Piilani Noguchi 05/12/2022

%Currently analyzes 2D data.
%Given cell image (and a b/w image showing where the cell area is) and orientation map, show:
% - local orientation varianece ie degree of parallel alignment between
%   between MTs
% - alignment with neatrest orientation ot nearest point on the cell border
%   (or, optionally, center)
% - histogram of of the alignment at the different distance ranges

%code commented out contained in the actual working program is removed here, so
%optional and alternative implemenatation is missing

%Orientation map is
%prepared  (currently using functions like BlkSVDOrient3D, all by Xiaoguang Feng)

% I'd like to imporove the terminology used in the comments to be
% consistent

%23-May-2018 added to input: a 2nd mask:mt image intensity threshold and 2nd image (insulin) to show its values in output according to new mask
%created a gui envelop with select buttons for 4 files and run

%%

%CLEARING VARIABLES & SETTING PATH


clearvars
% 
global diffsrgbBool0

if ispc
 path(path, 'Y:\Kaverina_Users\4136SDC_Piilani\MTorientation');
 path(path, 'Y:\Kaverina_Users\4136SDC_Piilani\MTorientation\CircStat2012a');
 path(path, 'Y:\Kaverina_Users\4136SDC_Piilani\MTorientation\Code\OMest');
 path(path,'Y:\Kaverina_Users\4136SDC_Piilani\MTorientation\Data');
else
path(path,'./CircStat2012a');
path(path,'./');
path(path,'./Data');
path(path,'./Code');
path(path,'./Code/OMest');
end
%%

%LOADING DATA 

%dataSet = 1; %which cell's photo, mask and orientation map files to use
constantNpixels = 5000;
noisePar=.0001;
calculate_local_angular_variance_swtich = false;
process_2nd_image = true;
if nargin ==0

%take input from UI
OMblocksize = 3; 
OMblockoverlap = 0;
nbins=4;%in degree of alignment histogram

%User grabbing files
         disp('Select image file');
         [Ifile, Idir]=uigetfile('*','Select image');      
         disp('Select mask file');
         [ maskfile ,maskdir]=uigetfile('*','Select mask');
   
        process_2nd_image = false;
        
elseif  nargin ~= 6
    
    error('Wrong number of input arguments');
    
else
   
    [Idir,Ifile,maskdir,maskfile,MTthreshold_dir,MTthreshold_file, I_2nd_dir, I_2nd_file] = ...
    deal(varargin{:});


%nbins=4;%in degree of alignment histogram



    %load new input data: MT image intensity threshold made in advance and 2nd image -
    % to match up against result orientation map
        I_2nd = imread([I_2nd_dir, I_2nd_file]);
        mask_MTthreshold = double(imread([ MTthreshold_dir, MTthreshold_file]));
   
end
      

 I = imread([Idir, Ifile]);
 mask = double(imread([ maskdir maskfile]));
   
        
% 

 if ~isa(I,'uint8')
     error('Input image must be 8 bit');%for now
     return;     
 end
 
 if ~all(size(I)==size(mask))
     
     error('Image and mask must me the same size.')
     return
     
 end
 
 
 %request mrightness threshold TholdToMaskAreaRatio, betwwen 0 and 1, 1
 %mening use all pixels as showing mocrotubules, 0 meanaing only the very brightest  (or none)
 
 %request micrometers per pixel value

micmperpixel = .03;%or .05
%if for some reazon input images are not grayscale
I=I(:,:,1);
mask=mask(:,:,1);




%%

%PREPARING DATA

fprintf(1, 'Loading orientation map... ');% here we presume image is ready - denoised etc.
I = double(I);

%OMname = [Ifile '_current_OM'];
%load('deconv_.tif_current_OM.mat');%for testing
OM=BlkSVDOrient(I,1,OMblocksize,OMblockoverlap);

fprintf(1, '  Done\n'); 


%May 2018: images anf masks now between 0 and one;
I_8bit = I;
I = I/255;
mask = mask/255;
%


%PrepDataSamples
%involves some manual approxiamtion

if process_2nd_image%process_2nd_image also meant a custom mt mask was given

    thold = mask_MTthreshold;
    npixesInData = sum(logical(mask_MTthreshold(:)).*logical(mask(:)));
    
else
%determine which pixels show an MT:
TholdToMaskAreaRatio =.3291;%custom threshold - which pixel, by brightness, show a mictotubule
npixesInMask = sum((mask(:)));
tholds=[0.01:.025:.5];
a1 = repmat(permute(tholds,([3 1 2])),size(mask));
a2 = repmat(I,[1,1,length(tholds)]) >a1;
repmask = repmat(mask,[1,1,length(tholds)])/255;%apply mask to a2
a3 = sum(sum(a2.*repmask));
a3=permute(a3,([1 3 2]));
a4= a3 ./npixesInMask - TholdToMaskAreaRatio;
[~,a5]=min(abs(a4));
thold = tholds(a5);%
npixesInData = sum(((I(:))>thold).*logical(mask(:)));


end



percentpixels = constantNpixels/npixesInData;
if percentpixels>1
    percentpixels = 1;
    disp('warining - too few pixels in data');
end




if ~(all(size(I)==size(OM)))
I=imresize(I,size(OM),'bilinear');%all loaded data images must be the same size
warning('Input images are not all the same size.');
end

mask = double(logical(mask));
if ~all(size(mask)==size(OM))
   double((imresize(mask,size(OM),'nearest'))) 
end
mask = 1-mask;%also, must me correct value range ()

%c=repmat(mask,[1,1,3]);

%%
%ROTATE DATA
if false%Below is one way of treating orientaion analysis. Currently not used, all prepared in advance, loaded above.
I0=I;
hists=[];
bins=[-pi/2 + (pi/360)/2 : pi/360: pi/2 - (pi/360)/2];
for Actr = 0:1:359

  I= imrotate(I0,Actr,'nearest','crop');
  OM=BlkSVDOrient(double(I),1,2,0);
  angs=atan2(real(OM),imag(OM));
  ex=(angs.*((I)>thold).*~mask);
  hists(end+1,:)=hist(ex(ex~=0),bins);
  fprintf(1, [num2str(Actr), ' ']);
end
figure
imagesc(hists);
return
end
fprintf(1, 'Data Rotated')
%%
%FIND MASK BOUNDARY

%determine pixels that make up the cell border:
[extmp ,wytmp] = meshgrid((1:size(I,2)),(1:size(I,1)));
maskinds1 =  wytmp(~logical(mask(:)));%coordinates of pixels that make up the cell
maskinds2 =  extmp(~logical(mask(:)));
maskindslinear = sub2ind(size(I),maskinds1,maskinds2);%1D indecies of above 2D coordinates 
boundscell=bwboundaries(mask)';%analysis to find borders in the mask image
cs=cell2mat(cellfun(@size,boundscell,'uni',false)');
[~,whichCell]=min(cs(:,1)); %bwboundaries result has two parts - actual cell border and a rectangular outline of the image
boundsinds = boundscell{2};                                         %so we want to chose the right one
boundsindsx = boundsinds(:,1);                                      %actually not guaranteed to work automatically
boundsindsy = boundsinds(:,2);                                      %check for each cell manually to be sure.
%because cell was outlained manually, the bounary is masy and crooked,
%smooth it:
smboundsx=smooth(boundsindsx,0.05,'loess');smboundsy=smooth(boundsindsy,0.05,'loess');
boundsO = atan(-diff(smboundsx)./diff(smboundsy))+pi/2;%orientation at points on border
boundsO(end+1)=boundsO(1)-boundsO(end);%add connecting point
boundsO=mod(boundsO,pi);%to wrap around connecting point;
%now create images from the above coordinate lists
boundsOmap = zeros(size(mask));
boundsmap = zeros(size(mask));
boundsindslinear=sub2ind(size(I),boundsindsx,boundsindsy);
boundsOmap(boundsindslinear)=boundsO;
boundsmap(boundsindslinear)=1;
bounDists = zeros(size(I));
fprintf(1, 'Mask Boundary found')
%%
%deal with orientation map, combine with mask, etc:
%OM comtains complex numbers, get the angle values from it
cm=colormap('hsv');
angs=atan2(real(OM),imag(OM));
angsNaN=angs;%remember to add to mask
angs2(isnan(angs))=0;
angs2=(angs-min(angs(:)));
angs2=(angs2/max(angs2(:)));
angs3=imresize(angs2,size(I),'nearest');
angsrgb=ind2rgb(floor(angs3.*63+1),cm);
I2 = I;%imresize(I,size(angs),'nearest');
I3 = repmat(I,[1 1 3]);
ex=((I3).*angsrgb);
ex2=rgb2ind(ex,64);
Omap = double(rgb2ind((I3).*angsrgb,64));%scrap
 
 %more set up including, neighborhood kernel
 nhoodKernel1=[ 1 2 3;1 2 3; 1 2 3];
 nhoodKernel2=[1 1 1;2 2 2;3 3 3];
 nhoodKernel12=nhoodKernel1-2;
 nhoodKernel22=nhoodKernel1-2;


 

  isnans=zeros(size(OM));
 Ithold = double(I2>thold);
 
 
 
 if calculate_local_angular_variance_swtich
     
     calculate_local_angular_variance();
     
end %end calculate_local_angular_variance

fprintf(1, 'OM done mask combined')


  %%
%FindDists
%find distances from nearest border point a given pixel 
distsFromC= zeros(size(I));
distsFromC0= zeros(size(I));
bounDists=zeros(size(I))-1;
nearestOs=zeros(size(I));
cellcenter= [mean(maskinds2) mean(maskinds1)];

fprintf(1, 'Calculating distance from border map... ');%... and the orientation of the neares point on boundary
 for indctr = 1:length(maskinds1)%for each pixel in cell
     
     distsTmp = [ boundsindsx boundsindsy] - repmat([maskinds1(indctr) maskinds2(indctr)],[length(boundsindsx),1]);                 
     [minVal, minI] = min(sum(distsTmp.^2,2).^.5);%actual distance
     bounDists(maskinds1(indctr), maskinds2(indctr))=minVal;
     
     %we also look at oter distances: from cell center
     distsFromC0(maskinds1(indctr), maskinds2(indctr)) =...
         ((maskinds2(indctr)-cellcenter(1)).^2 + (maskinds1(indctr)-cellcenter(2)).^2).^.5;
     
     distsFromC(maskinds1(indctr), maskinds2(indctr)) =...
         ((boundsindsy(minI)-cellcenter(1)).^2 + (boundsindsx(minI)-cellcenter(2)).^2).^.5;
     
     %orientation at the nearest boundary pixel:
     nearestOs(maskindslinear(indctr)) = boundsOmap(boundsindslinear(minI));
 end
fprintf(1, '  Done\n');
 
 
fprintf(1, 'Binning and analysis... ');
%finalize the result
angs= (angs+pi/2);
angs= mod(angs+pi/2,pi);
 nearestOs = mod(nearestOs+pi/2,pi);
diffs=angs-nearestOs;
 

draw_masked_angle_delta_image();%create good image
 
 
%% from here on various outputs wwill happen
Ifile2 = Ifile(1:end-5);%cut off the extension from filename to create output file name
filename_timestamp = datevec(now);%all output files will have a timestamp unique for each run
filename_timestamp(end)=floor(filename_timestamp(end));
filename_timestamp=filename_timestamp(2:end);
filename_timestamp=num2str(filename_timestamp);
filename_timestamp(filename_timestamp==' ') = '';

%% Show the map of local angles relative to that of border
%parallel to border is one extreme, perpendicular the other
%create a custom colormap to show background in gray, data pixels in color
%and, 2018, boundary in black
ex2=((diffs).*~mask.*Ithold);
ex2(ex2>pi/2)= ex2(ex2>pi/2)-pi;
ex2(ex2<-pi/2)= ex2(ex2<-pi/2)+pi;
ex2=abs(ex2);
diffsReady= ex2;
cmj=colormap('jet');
cmj(end,:)=[0 0 0]; %background color
cmj(1,:)=1; %separate color for border% todo: ensure unique color ie white
ex2(ex2==0)=max(ex2(:))+.1;

ex2((logical(boundsOmap))) =  0;
deltaangleFigure = figure;

outputimage_ex2 = floor(ex2/(.5*pi)*64);
deltaangleFigure_image = imagesc(outputimage_ex2);
colormap(cmj);

%write image with custom colormap to file1
exrgb = ind2rgb(outputimage_ex2,cmj);
filename = [Ifile2 '_imageout' filename_timestamp    '.tif'];

imwrite(exrgb,filename);

%savefigure %image just in case - currently to use colormap


%%
%FindDataBins
%Find degrees of alignment for <nbins> areas of different offsets form border ranges

%Initiate variables.
binSubject=diffsReady;
bounDists(isinf(bounDists))=0;
circumbin=[];
circumbinNs=[];
maxdist=max(bounDists(:));
circboundimg=0;


binarea = sum((diffsrgbBool0(:)))/nbins;% histogram bins each contain <binarea> values

%this was in previous version... why?
%binarea = sum((diffsrgbBool0(:)))/nbins;% histogram bins each contain <binarea> values

ex1=binSubject(diffsrgbBool0);
ex2=bounDists(diffsrgbBool0);

ex1_degrees = 180 * ex1/(pi);

[exhistout1, exhistout2] = hist3([ex2,ex1_degrees],[4,10]);%flipped to y, x!!!!!!!!!!!!!
%[exhistout1, exhistout2] = hist3([ones(length(ex1),1),ex2]);%[angledelta , distance]


%norm across distance bins
exhistout1_normd = exhistout1 ./ repmat(sum(exhistout1,2),[1,10]);
exhistout1_normd = exhistout1_normd*100;%percent

figure;
imgtmp = imagesc(exhistout1_normd);
figure;
imgtmp = imagesc(exhistout1);



% 
% set(imgtmp,'XData',exhistout2{2})
% set(imgtmp,'YData',exhistout2{1})

set(gca,'YDir','normal');

%..... 
h = colorbar;

set(get(h,'title'),'string','...','Rotation',90.0);

filename=[Ifile2 '_Dsisp3Dout'];
filename=[filename filename_timestamp];

[ex2_sorted, ex2_sorted_inds]=sort(ex2);
ex1_sorted=ex1(ex2_sorted_inds);
% 
% ex1_sorted_str = ['angle, redians '  num2str(ex1_sorted')];
% ex2_sorted_str = ['distance, pixels '  num2str(ex2_sorted')];%

extable=table(num2str(ex1_sorted), num2str(ex2_sorted));
extable.Properties.VariableNames={'angle' 'distance'};
writetable(extable,[filename '.txt']);
disp([filename ' saved'])






%return
%begin numeric analysis

rctr = 1;
hibound=1;lowbound=0; %values are number of pixels, initially range is between 0 and one pixel from border
binN=1;

%In order to fill the <nbins> bins neatly and equally, we iterate gradually through
%groups of pixels located wthin a range of  <hibound> to <hibound> pixels
%from border (or optionally center). Is this overkill?


%direct comparison: ex=(diffsrgbBool0.*bounDists);...;[exout1 exout2] = hist3([ex1,ex2n])

while hibound<maxdist

    hibound=hibound+.1; %<hibound> is raised gradually untill the number of pixels in range
                        %matches the bin size (<binarea>).
                        %Then bin is created and <lowbound> is mover to <hibound>
    circumbintmp =bounDists<hibound&bounDists>lowbound;
    circumbintmp = circumbintmp&diffsrgbBool0;
    binareatmp =  sum(circumbintmp(:));
    if binareatmp>=binarea
        lowbound=hibound;        
        
        %Take only <percentpixels> percent of pixels from zone.
        %This is a step towards normalizing results between numerous cells,
        %which have different areas.
        bintmp = binSubject(circumbintmp(:)&diffsrgbBool0(:));
        binsampleinds = randperm(length(bintmp));
        binsampletmp=bintmp(binsampleinds(1: floor(length(bintmp)*percentpixels)));
        circumbin =[circumbin; binsampletmp];
        
        circumbinNs = [circumbinNs; ones(length(binsampletmp),1)*binN ];
        binN=binN+1;
%        circboundimg = circboundimg+ circumbintmp*hibound.*binSubject.*diffsrgbBool0;
        circboundimg = circboundimg+ circumbintmp*hibound.*binSubject.*diffsrgbBool0;
        
        bin_boundary=bounDists<(lowbound) & bounDists>(lowbound-3);
        %show bins in circboundimg
        
        circboundimg(bin_boundary) = 0;
    
    end   
end
%make the final, INNERMOST (#<nbins>) bin  (special case)
        bintmp = binSubject(circumbintmp(:)&diffsrgbBool0(:));
        binsampleinds = randperm(length(bintmp));
        binsampletmp=bintmp(binsampleinds(1: floor(length(bintmp)*percentpixels)));
        circumbin =[circumbin; binsampletmp];
        circumbinNs = [circumbinNs; ones(length(binsampletmp),1)*binN ];
 circboundimg = circboundimg+ circumbintmp*hibound.*binSubject.*diffsrgbBool0;
    
 
 %%
%AnalyzeAndSave
[P,ANOVATAB,STATS] =anova1(circumbin/(pi/2),(circumbinNs));
figure
[c,m,h,nms] = multcompare(STATS);

datasetOut = struct('mask',mask,'I',I,'binSubject',binSubject,'angs',angs,'circumbin',circumbin,'circumbinNs',circumbinNs,'multOut',{c,m,h,nms});
%num2str(dataSet)
datasetName=[Ifile '_DatasetOutRadDiff2' num2str(dataSet) '.mat'];
disp(datasetName);
disp(size(circumbinNs));
save(datasetName,'datasetOut');

return







 function calculate_local_angular_variance() 
        nhoodsize = 8;
        
 ItholdPadded = padarray(Ithold,[nhoodsize nhoodsize],'circular');
 Ipadded = padarray(I2,[nhoodsize nhoodsize],'circular');
        
 tmpIthold = zeros(size(OM));
 
  OmapPadded = padarray(2*angs,[nhoodsize nhoodsize],'circular');
 OmapPadded(isnan(OmapPadded))=0;
 
 [x,y]=meshgrid(-nhoodsize:nhoodsize,-nhoodsize:nhoodsize);
 circKernel = sqrt(x.^2+y.^2)<nhoodsize;
% circKernel(nhoodsize:nhoodsize+2,nhoodsize:nhoodsize+2)=0;
tmpItholdCum= zeros(size(OM));
divisorCum=zeros(size(OM));
 angVarmat = zeros(size(OM));
%%
fprintf(1, 'Calculating local angular variance... ');
%angular variance is simply the sum of angles, larger means more of
%a particular direction, smaller - more evenly distributed directions
  for ctr1=-nhoodsize:nhoodsize
      for ctr2=-nhoodsize:nhoodsize
          
          %not part of clean comments, but better leave this comment here, might be important:
          % threshold values!!!!!!!!!!!!!!!!!!!!!!!!
        if circKernel(ctr1+nhoodsize+1,ctr2+nhoodsize+1)%shape of neighborhood may be square or circular, the latter means we skip some pixels
      %see if give neighbor pixel is part of an MT, use only those which
      %are
      tmpIthold = ItholdPadded((nhoodsize+1:end-nhoodsize)+ctr1,(nhoodsize+1:end-nhoodsize)+ctr2);
      tmpItholdCum = tmpItholdCum+tmpIthold;
      divisortmp = Ipadded((nhoodsize+1:end-nhoodsize)+ctr1,(nhoodsize+1:end-nhoodsize)+ctr2);
      divisorCum = divisorCum + divisortmp;
      tmp=OmapPadded((nhoodsize+1:end-nhoodsize)+ctr1,(nhoodsize+1:end-nhoodsize)+ctr2);
   
      angVarmattmp =  tmpIthold .* exp(1i.*tmp);
      diffTmp = angs-tmp;
      morethanpi =  diffTmp>pi/2;
      lessthanmpi = diffTmp < -pi/2;
      diffTmp(morethanpi) =diffTmp(morethanpi)- pi;
      diffTmp(lessthanmpi) =diffTmp(lessthanmpi)+ pi;
      angVarmat=angVarmat+angVarmattmp;
        end
      end
  end
  
  %now we're curious about the result, so check out the result
  %might want to look at more than one version of it
  ex=(abs(angVarmat(:))./(tmpItholdCum(:)+1));
  meanAngVarmat= mean(ex(logical(Ithold)&~mask));
  fprintf(1, '  Done\n');
  tmpIthold3=repmat(tmpItholdCum,[1,1,3]);
  angVarmat0=angVarmat;
  angVarmat=abs(angVarmat)./(tmpItholdCum+1);
  %figure;imagesc(((I3).*angsrgb).*repmat(1-mask,[1,1,3]))
  
  
        
    end





    function draw_masked_angle_delta_image()
       
%create masked diff_angle color map
cm2=cm;
cm2(1,:)=[0 0 0];%create colormap for masked angle map from hsv by making first color black 
diffsrgb=ind2rgb(floor((abs(diffs))/(pi).*63+1),cm2);
%figure;
%%count # of pixels in data, those in diffsrgbBool0==1
%old syntax - needs purging
diffsrgb3=(((I3>thold).*diffsrgb).*repmat(1-mask,[1,1,3]));
diffsrgbBool=(sum(diffsrgb3,3)==0);
diffsrgbBool0=~logical(sum(diffsrgbBool,3));
diffsrgbBool=repmat(diffsrgbBool,[1,1,3]);
diffsrgb3(diffsrgbBool)=1;
 imagesc(diffsrgb3);
% colormap(cm2);
masked_diffsrgb3 = ((I3>.3).*diffsrgb).*repmat(1-mask/255,[1,1,3]);
 imagesc(masked_diffsrgb3);%good image
%end create masked diff_angle color map 
    end


end