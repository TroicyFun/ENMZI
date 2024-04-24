clear

run def_rg_enmzi.m
mth0 = 5:9;
LM0 = length(mth0);
mth = 6:8;
[~,iam,~] = intersect(mth0,mth);
nrg = size(rg_enmzi,1);
end_year = 2023;

recname = '东亚季风北边缘带指数（ENMZI）.xlsx';
HF = size(dir(recname));
if HF(1)
    delete(recname)
end

iy = 1;
for I = 1:4
if I==1  % ERA5
dsnm = 'ERA5';
yr = 1940:end_year;
LY = length(yr);
enmzi0 = nan(LY,LM0+1,nrg);
vrnm = 'uwnd';
vrnm_rd = 'u';
for i = 1:LY
    fn = ['F:\ERA\era5\monthly\',vrnm,'\',vrnm,'_',num2str(yr(i)),'.nc'];
    if i==1
        lon = ncread(fn,'longitude');
        lat = ncread(fn,'latitude');
        lev = ncread(fn,'level');
    end
    for a = 1:nrg
        [ ISR ] = get_isr( lon,lat,rg_enmzi(a,:,1),1 );
        [LON,LAT] = meshgrid(lon(ISR(1):ISR(2)),lat(ISR(4):ISR(3))); 
        lv0 = lv_enmzi;
        [~,IS0] = min(abs(lev-lv0));
        if i~=LY
            tc = squeeze(ncread(fn,vrnm_rd,[ISR(1) min(ISR(3:4)) IS0 mth0(1)],[ISR(2)-ISR(1)+1 abs(diff(ISR(3:4)))+1 1 LM0]));
        else
            tc = squeeze(ncread(fn,vrnm_rd,[ISR(1) min(ISR(3:4)) IS0 1 mth0(1)],[ISR(2)-ISR(1)+1 abs(diff(ISR(3:4)))+1 1 1 LM0]));
        end
        tc(:,:,LM0+1) = mean(tc(:,:,iam),3);
        for m = 1:LM0+1
            tc_c = tc(:,:,m)';
            LAT_c = LAT;
            LAT_c(isnan(tc_c)) = nan;
            enmzi0(i,m,a) = squeeze(sum(sum(tc_c.*cosd(LAT_c),1,'omitnan'),2,'omitnan')./sum(sum(cosd(LAT_c),1,'omitnan'),2,'omitnan'));
        end
    end
end

elseif I==2  % JRA55
dsnm = 'JRA-55';
yr = 1958:end_year;
LY = length(yr);
vrnm = 'ugrd';
lv0 = lv_enmzi;
enmzi0 = nan(LY,LM0+1,nrg);
for i = 1:LY
    yr_nm = num2str(yr(i));
    fn = ['F:\RD\JRA-55\monthly\plevels\',vrnm,'\anl_p125.033_',vrnm,'.',yr_nm,'01_',yr_nm,'12'];
    nc = ncgeodataset(fn);
    nv = nc.variables;
    if i==1
        lon = nc.data('lon');
        lat = nc.data('lat');
        lev = nc.data('isobaric');
        [~,IS] = min(abs(lev-lv0));
    end
    t0 = squeeze(nc.data(nv{2}));
    for a = 1:nrg
        [ ISR ] = get_isr( lon,lat,rg_enmzi(a,:,1),1 );
        [LON,LAT] = meshgrid(lon(ISR(1):ISR(2)),lat(ISR(4):ISR(3))); 
        tc = permute(squeeze(t0(mth0,IS,ISR(4):ISR(3),ISR(1):ISR(2))),[2 3 1]);
        tc(:,:,LM0+1) = mean(tc(:,:,iam,:),3);
        for m = 1:LM0+1
            tc_c = tc(:,:,m);
            LAT_c = LAT;
            LAT_c(isnan(tc_c)) = nan;
            enmzi0(i,m,a) = squeeze(sum(sum(tc_c.*cosd(LAT_c),1,'omitnan'),2,'omitnan')./sum(sum(cosd(LAT_c),1,'omitnan'),2,'omitnan'));
        end
    end
end

elseif I==3 || I==4 % NCEP/NCAR, 20CR
vrnm = 'uwnd';
if I==3
    dsnm = 'NCEP/NCAR';
    yr = 1948:end_year;
    fn = ['F:\RD\NCEPNCAR\',vrnm,'.mon.mean.nc'];
elseif I==4
    dsnm = '20CR';
    yr = 1836:2015;
    fn = ['F:\RD\20CR\',vrnm,'.mon.mean.nc'];
end
LY = length(yr);
enmzi0 = nan(LY,LM0+1,nrg);
lon = double(ncread(fn,'lon'));
lat = double(ncread(fn,'lat'));
lev = double(ncread(fn,'level'));
ik = 0;
for i = yr
    tn(ik+1:ik+LM0) = (i-yr(1))*12+mth0;
    ik = ik+LM0;
end
for a = 1:nrg
    [ ISR ] = get_isr( lon,lat,rg_enmzi(a,:,1),1 );
    [LON,LAT] = meshgrid(lon(ISR(1):ISR(2)),lat(min(ISR(3:4)):max(ISR(3:4)))); 
    lv0 = lv_enmzi;
    [~,IS0] = min(abs(lev-lv0));
    t0 = squeeze(ncread(fn,vrnm,[ISR(1) min(ISR(3:4)) IS0 1],[ISR(2)-ISR(1)+1 abs(diff(ISR(3:4)))+1 1 inf]));
    [nx,ny,~] = size(t0);
    tc = reshape(t0(:,:,tn),nx,ny,LM0,LY);
    tc(:,:,LM0+1,:) = mean(tc(:,:,iam,:),3);
    for i = 1:LY
        for m = 1:LM0+1
            tc_c = tc(:,:,m,i)';
            LAT_c = LAT;
            LAT_c(isnan(tc_c)) = nan;
            enmzi0(i,m,a) = squeeze(sum(sum(tc_c.*cosd(LAT_c),1,'omitnan'),2,'omitnan')./sum(sum(cosd(LAT_c),1,'omitnan'),2,'omitnan'));
        end
    end
end
end

enmzi = squeeze(sum(permute(repmat(sg_enmzi,LY,1,LM0+1),[1 3 2]).*enmzi0,3));
idn = ['ENMZI - ',dsnm];
wt_cal = {idn,'M','C5C9','m/s'};
xlswrite(recname,wt_cal,1,['A',num2str(iy)]);
xlswrite(recname,{[num2str(LY),'\m']},1,['A',num2str(iy+1)]);
xlswrite(recname,1:12,1,['B',num2str(iy+1)]);
xlswrite(recname,{'JJA-mean'},1,[char('A'+12+1),num2str(iy+1)]);
xlswrite(recname,yr',1,['A',num2str(iy+2)]);
xlswrite(recname,enmzi(:,1:5),1,[char('A'+mth0(1)),num2str(iy+2)]);
xlswrite(recname,enmzi(:,6),1,[char('A'+12+1),num2str(iy+2)]);
iy = iy+LY+3;
end

