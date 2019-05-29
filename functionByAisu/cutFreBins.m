function [freIndex] = cutFreBins(f,freRange) 

	idxStart = find(f>freRange(1),1,'first');
	idxEnd   = find(f<freRange(2),1,'last');
	freIndex = [idxStart:idxEnd];