function diagMat = doDiags(aMat)
%Given Matrix fullMat, gets diags for t_years
    size_ = size(aMat, 1);
    t_years = 150 + size_ - 1;
    diagMat = zeros(size_, t_years);
    
    temp = [repmat(aMat(:, 1), 1, size_ - 1), aMat, repmat(aMat(:, end), 1, size_ - 1)];

    start_ = size_ - 1;
    for i = 1:t_years
        ind1 = i;
        ind2 = start_ + i;
        diagMat(:, i) = diag(temp(:, ind1:ind2));
    end

    if t_years < 229
        diagMat = [repmat(diagMat(:, 1), 1, 79) diagMat repmat(diagMat(:, end), 1, (229-t_years))];
    else
        diagMat = [repmat(diagMat(:, 1), 1, 79) diagMat];
    end

    