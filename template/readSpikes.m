function [tspI,l] = readSpikes(DIR,n,filename)
    filepath = [DIR,filename];
    f1 = fopen(filepath);
    tspI = struct([]);
    data = fread(f1,[2,inf],'float');
    fclose(f1);
	if ~ sum(size(data))==0
		ispike = data(1,:);
		tspike = data(2,:);
		for i=1:n
			tspI(i).tsp = tspike(ispike==i);
			l(i) = sum(ispike==i);
		end
	else
		tspI(n).tsp = [];
		l(n) = 0;
	end
end
