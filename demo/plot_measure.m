function plot_measure(filename,measure)
    if exist(filename, 'file') ~= 2
        fprintf("File does not exist.");
        return;
    end
  
    coor = extract_coor(filename);
    if strcmp(measure,'dirac')  || strcmp(measure,'d')
        plot(coor(:,1),coor(:,2),'*','linewidth',5,'Markersize',5);
    elseif strcmp(measure,'empirical')  || strcmp(measure,'e')
        plot(coor(:,1),coor(:,2),'.');
    end
       
end

function coor = extract_coor(filename)
    fid = fopen(filename,'r');
    txt = textscan(fid,'%s','delimiter','\n'); 
    fclose(fid);
    
    txt = txt{1};
    header = strsplit(string(txt{1}),' ');
    dim = str2double(header(3));
    num = str2double(header(2));
    
    coor = zeros(num,dim);
    
    for i = 2:num+1
       tmp = extractAfter(string(txt{i}),"}");
       coor(i-1,:)=str2double(regexp(tmp,'[+-]?\d?\.?\d+','match')');
    end
end
