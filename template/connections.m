function [preStruct,preE,preI] = connections(theme,lgnfile,ntotal,excite)
    load(lgnfile);
    theta = [reshape(etheta,[p.nv1e,1]); reshape(itheta,[p.nv1i,1])];
    fid = fopen([theme,'/','cMatrix_summary']);
    postE = fread(fid, [1,ntotal],'int32');
    postI = fread(fid, [1,ntotal],'int32');
    preE = fread(fid, [1,ntotal],'int32');
    preI = fread(fid, [1,ntotal],'int32');
    mat = fread(fid, [ntotal,ntotal],'int8');
    preStruct = struct([]);
    for i=1:ntotal
        preStruct(i).preE = preE(i);
        if preE(i)>0
            id = find(mat(i,:) >= 1 & excite > 0.5, preE(i),'first'); 
            preStruct(i).ID_exc = id;
            preStruct(i).s_exc = mat(i,id);
            preStruct(i).theta_exc = theta(preStruct(i).ID_exc);
            preStruct(i).nLGN_exc = nLGN(preStruct(i).ID_exc,1);
        else 
            preStruct(i).ID_exc = [];
            preStruct(i).theta_exc = [];
            preStruct(i).nLGN_exc = [];
        end
        preStruct(i).preI = preI(i);
        if preI(i)>0
            preStruct(i).ID_inh = find(mat(i,:) >= 1 & excite < 0.5, preI(i),'first'); 
            preStruct(i).theta_inh = theta(preStruct(i).ID_inh);
            preStruct(i).nLGN_inh = nLGN(preStruct(i).ID_inh,1);
        else
            preStruct(i).ID_inh = [];
            preStruct(i).theta_inh = [];
            preStruct(i).nLGN_inh = [];
        end
        preStruct(i).LGN_onid = v1Map{i,1};
        preStruct(i).LGN_offid = v1Map{i,2};
    end
    fclose(fid);
end
