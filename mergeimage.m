function [Intensity_merged, RGB_merged, FRET_merged]=mergeimage(D_image, A_image, map)

[X_size, Y_size]=size(D_image);
Intensity_merged=zeros(X_size,Y_size);
FRET_merged=zeros(X_size,Y_size);

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
            if mappedx < 1, mappedx = 1; end
            if mappedy < 1, mappedy = 1; end
            if mappedx > X_size, mappedx = X_size; end
            if mappedy > Y_size, mappedy = Y_size; end
            Intensity_merged(x,y)=(D_image(x,y)+A_image(mappedx,mappedy))/2;
            FRET_merged(x,y)=A_image(mappedx,mappedy)/(D_image(x,y)+A_image(mappedx,mappedy));
            mapped_red(x,y)=A_image(mappedx,mappedy);
        end
    end
else
    for x=1:X_size
        for y=1:Y_size
            [mappedx, mappedy] = map_pt2pt(x,y,map);
%             mappedx= round(map(1,1) +...
%                 map(1,2)*x + map(1,3)*y +...
%                 map(1,4)*x^2 + map(1,5)*y^2 +...
%                 map(1,6)*x^3 + map(1,7)*y^3 +...
%                 map(1,8)*x*y +...
%                 map(1,9)*x^2*y + map(1,10)*x*y^2);
%             mappedy= round(map(2,1) +...
%                 map(2,2)*x + map(2,3)*y +...
%                 map(2,4)*x^2 + map(2,5)*y^2 +...
%                 map(2,6)*x^3 + map(2,7)*y^3 +...
%                 map(2,8)*x*y +...
%                 map(2,9)*x^2*y + map(2,10)*x*y^2);

            mappedy=mappedy-Y_size;
            if mappedx < 1, mappedx = 1; end
            if mappedy < 1, mappedy = 1; end
            if mappedx > X_size, mappedx = X_size; end
            if mappedy > Y_size, mappedy = Y_size; end
            Intensity_merged(x,y)=(D_image(x,y)+A_image(mappedx,mappedy))/2;
            FRET_merged(x,y)=A_image(mappedx,mappedy)/(D_image(x,y)+A_image(mappedx,mappedy));
            mapped_red(x,y)=A_image(mappedx,mappedy);
        end
    end
end
bg=median(mapped_red(:));
RGB_merged(:,:,1)=(mapped_red-bg)/(max(mapped_red(:))-bg);
bg=median(D_image(:));
RGB_merged(:,:,2)=(D_image-bg)/(max(D_image(:))-bg);
RGB_merged(:,:,3)=zeros(X_size, Y_size);

