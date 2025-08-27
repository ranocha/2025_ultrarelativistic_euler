function col = grb(x)

	if (0 <= x) && (x <= 0.2)

		z = 4*x;
		a = 1;
		b = 1;
		c = 1-z;

	elseif (0.2 <= x) && (x <= 0.3)

		z = 10*x -2;
		a = 1;
		b = -0.2*z*z +1;
		c = 0.2*z*z -0.4*z +0.2;

	elseif (0.3 <= x) && (x <= 0.45)

		z = 4*x -1;
		a = 1;
		b = 1-z;
		c = 0;

	elseif (0.45 <=x) && (x <= 0.55)

		z = 10*x -4.5;
		a = -0.2*z*z +1;
		b = 0.2*z*z -0.4*z +0.2;
		c = 0.2*z*z;

	elseif (0.55 <= x) && (x <= 0.7)

		z = 4*x -2;
		a = 1-z;
		b = 0;
		c = z;

	elseif (0.7 <= x) && (x <= 0.8)

		z = 10*x -7;
		a = 0.2*z*z -0.4*z +0.2;
		b = 0;
		c = -0.4*z*z +0.4*z + 0.8;

	elseif (0.8 <= x) && (x <= 1)

		z = 4*x-3;
		a = 0;
		b = 0;
		c = 1-z;

	end

	col = [a,b,c];
end
