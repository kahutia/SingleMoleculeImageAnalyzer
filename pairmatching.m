function [pairsfound, num_paired]=pairmatching(mol_pos, num_mol, map)
num_paired=0;
pairsfound=zeros(num_mol,4);
pos_uncertainty=3;

for i=1:num_mol
    if mol_pos(i,2) < 250
        [~,N_order]=size(map);
        N_order=(N_order-1)/2;
        
        if N_order==3
            mappedx= round(map(1,1) +...
                map(1,2)*mol_pos(i,1) + map(1,3)*mol_pos(i,2) +...
                map(1,4)*mol_pos(i,1)^2 + map(1,5)*mol_pos(i,2)^2 +...
                map(1,6)*mol_pos(i,1)^3 + map(1,7)*mol_pos(i,2)^3);
            mappedy= round(map(2,1) +...
                map(2,2)*mol_pos(i,1) + map(2,3)*mol_pos(i,2) +...
                map(2,4)*mol_pos(i,1)^2 + map(2,5)*mol_pos(i,2)^2 +...
                map(2,6)*mol_pos(i,1)^3 + map(2,7)*mol_pos(i,2)^3);
        else
            mappedx= round(map(1,1) +...
                map(1,2)*mol_pos(i,1) + map(1,3)*mol_pos(i,2) +...
                map(1,4)*mol_pos(i,1)^2 + map(1,5)*mol_pos(i,2)^2 +...
                map(1,6)*mol_pos(i,1)^3 + map(1,7)*mol_pos(i,2)^3 +...
                map(1,8)*mol_pos(i,1)^4 + map(1,9)*mol_pos(i,2)^4);
            mappedy= round(map(2,1) +...
                map(2,2)*mol_pos(i,1) + map(2,3)*mol_pos(i,2) +...
                map(2,4)*mol_pos(i,1)^2 + map(2,5)*mol_pos(i,2)^2 +...
                map(2,6)*mol_pos(i,1)^3 + map(2,7)*mol_pos(i,2)^3 +...
                map(2,8)*mol_pos(i,1)^4 + map(2,9)*mol_pos(i,2)^4);
        end
        % check if there is acceptor at the mapped position
        for k=1:num_mol
            if (mol_pos(k,1) > mappedx-pos_uncertainty && ...
                    mol_pos(k,1) < mappedx+pos_uncertainty && ...
                    mol_pos(k,2) > mappedy-pos_uncertainty && ...
                    mol_pos(k,2) < mappedy+pos_uncertainty )
                num_paired=num_paired+1;
                pairsfound(num_paired,1)=mol_pos(i,1);
                pairsfound(num_paired,2)=mol_pos(i,2);
                pairsfound(num_paired,3)=mol_pos(k,1);
                pairsfound(num_paired,4)=mol_pos(k,2);
            end
        end
    end
end

pairsfound=pairsfound(1:num_paired,:);
end