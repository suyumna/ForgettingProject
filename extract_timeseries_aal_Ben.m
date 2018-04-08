function time_series = extract_timeseries_aal_Ben(data, template)

% extracts the average time courses of BOLD signal from each of the regions
% defined in the AAL Tzourio-Mazoyer atlas. Data  should be 4D

dim = size(data);
dim(4)=1; % --> I put this thing here not to modify the previous script..silly I know --> our data are 3D
data = reshape( data, dim(1)*dim(2)*dim(3),dim(4) ); % [Voxels X Time] 
template = reshape( template, dim(1)*dim(2)*dim(3),1 ); % [Voxels X 1]
%  data2= reshape( data,[],1 );  --> Ben
%template = reshape( template,[],1 ); --> Ben


for i = 1:max(max(max(template)));
    x= find( template == i);
    temp=zeros(1,dim(4));
    for j = 1:length(x)
       % size( temp)
       % size( data(x(j),:))
        temp=temp+data(x(j),:);
    end
    time_series(i,:)=temp./length(x);
end

% time_series [ROI X Time]
% R=corr(time_series); % [ROI X ROI]
    
%  ============== THis is the original script

% function time_series = extract_timeseries_aal(data, template)
% 
% % extracts the average time courses of BOLD signal from each of the regions
% % defined in the AAL Tzourio-Mazoyer atlas.
% 
% dim = size(data);
% 
% data = reshape( data, dim(1)*dim(2)*dim(3),dim(4) );
% template = reshape( template, dim(1)*dim(2)*dim(3),1 );
% 
% for i = 1:max(max(max(template)));
%     x= find( template == i);
%     temp=zeros(1,dim(4));
%     for j = 1:length(x)
%        % size( temp)
%        % size( data(x(j),:))
%         temp=temp+data(x(j),:);
%     end
%     time_series(i,:)=temp./length(x);
% end
