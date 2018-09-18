function np_array = m2py(M, np)
	m = size(M,1);
	n = size(M,2);
	np_array = np.zeros({int32(m), int32(n)});
	for i = 1:m
		for j = 1:n
			np_array.itemset(int32(i-1), int32(j-1), double(M(i,j)));
		end
	end
end
