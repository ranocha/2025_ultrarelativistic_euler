function colmap = grb_colormap(N)

	colmap = zeros(N,3);

	for i = 1:N
	colmap(i,:) = grb((i-1)/(N-1));
	end

end
