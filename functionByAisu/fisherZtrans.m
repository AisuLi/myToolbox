function  [z] = fisherZtrans(r1,r2,n1,n2)





	for i = 1:length(r1) 
	    zf1(i) = 1/2*ln((1+r1(i))/(1-r1(i)));
		zf2(i) = 1/2*ln((1+r2(i))/(1-r2(i)));
		z(i) = (zf1(i)-zf2(i))/sqrt(1/(n1-3)+1/(n1-3));
	end
	


	meanR1 = mean(r1);
	meanR2 = mean(r2);

    zf1(i) = 1/2*ln((1+r1(i))/(1-r1(i)));
	zf2(i) = 1/2*ln((1+r2(i))/(1-r2(i)));
	z = (zf1(i)-zf2(i))/sqrt(1/(n1-3)+1/(n1-3));

	
