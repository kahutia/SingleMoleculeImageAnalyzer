function [D2Amap,A2Dmap]=calc_map_matrix(X_size, Y_size,map)
%% loading map coefficient from a map file
% [mapfilename mappathname]=uigetfile('*.mapcoef', 'mapping file');
% mapfid=fopen([mappathname mapfilename],'r');
% map=zeros(2,7);
% for j=1:2
%     for i=1:7
%         fscanf(mapfid,'%s',1);
%         map(j,i)=fscanf(mapfid,'%g',1);
%     end
% end
% fclose(mapfid);

if isempty(map)
    D2Amap=[];
    A2Dmap=[];
else
    %% Create Donor to Acceptor full map matrix
    D2Amap=zeros(X_size,Y_size,2);

    [~,N_order]=size(map);
    N_order=(N_order-1)/2;
    
    if N_order==3
        for x=1:X_size
            for y=1:Y_size
                mappedx= round(map(1,1) +...
                    map(1,2)*x + map(1,3)*y +...
                    map(1,4)*x^2 + map(1,5)*y^2 +...
                    map(1,6)*x^3 + map(1,7)*y^3);
                mappedy= round(map(2,1) +...
                    map(2,2)*x + map(2,3)*y +...
                    map(2,4)*x^2 + map(2,5)*y^2 +...
                    map(2,6)*x^3 + map(2,7)*y^3);
                mappedy=mappedy-Y_size;

                if mappedx < 1 || mappedx > X_size
                    mappedx=0;
                end
                if mappedy < 1 || mappedy > X_size
                    mappedy=0;
                end

                D2Amap(x,y,1)=mappedx;
                D2Amap(x,y,2)=mappedy;
            end
        end
    else
        for x=1:X_size
            for y=1:Y_size
                [mappedx, mappedy] = map_pt2pt(x,y,map);
%                 mappedx= round(map(1,1) +...
%                     map(1,2)*x + map(1,3)*y +...
%                     map(1,4)*x^2 + map(1,5)*y^2 +...
%                     map(1,6)*x^3 + map(1,7)*y^3 +...
%                     map(1,8)*x^4 + map(1,9)*y^4);
%                 mappedy= round(map(2,1) +...
%                     map(2,2)*x + map(2,3)*y +...
%                     map(2,4)*x^2 + map(2,5)*y^2 +...
%                     map(2,6)*x^3 + map(2,7)*y^3 +...
%                     map(2,8)*x^4 + map(2,9)*y^4);
                mappedy=mappedy-Y_size;
                
                if mappedx < 1 || mappedx > X_size
                    mappedx=0;
                end
                if mappedy < 1 || mappedy > X_size
                    mappedy=0;
                end
                
                D2Amap(x,y,1)=mappedx;
                D2Amap(x,y,2)=mappedy;
            end
        end
    end
    %% Create Acceptor to Donor full map matrix
    A2Dmap=zeros(X_size, Y_size);

    for x=1:X_size
        for y=1:Y_size
            if D2Amap(x,y,1) > 0 && D2Amap(x,y,1) <= X_size && ...
                    D2Amap(x,y,2) > 0 && D2Amap(x,y,2) <= Y_size 
                A2Dmap(D2Amap(x,y,1),D2Amap(x,y,2),1)=x;
                A2Dmap(D2Amap(x,y,1),D2Amap(x,y,2),2)=y;
            end
        end
    end
end