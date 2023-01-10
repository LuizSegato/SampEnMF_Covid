%This function calculates the SampEnMF curves for a given image

% [5] L. F. S. dos Santos, L. A. Neves, G. B. Rozendo, M. G. Ribeiro, M. Z. do Nascimento, T. A. A.
%Tosta, Multidimensional and fuzzy sample entropy (sampenmf) for quantifying h&e histological
%images of colorectal cancer, Computers in biology and medicine 103 (2018) 148–160.

% Input:
% features - SampEnMF values (attributes) provided from the SampEnMF function (first composition of vectors according [5])

% Output: 
% results - a feature vector composed by attributes A, Γ, O, W and E (second composition of vectors according [5]) 

function [ results ] = curve_metrics( features, k_start, k_increment, k_limit )
    
    x = k_start:k_increment:k_limit;
      
    %Eliminate unrevealing scales (with Inf e NaN)
    infs = isinf(features);
    nans = isnan(features);
    
    indices = ones(1, size(features,2));
    
    for j = 1 : size(features,2)
        for i = 1 : size(features,1)
            if infs(i, j) == true || nans(i, j) == true
                indices(1, j) = 0;
            end
        end
    end
    
    indices = logical(indices);

    features = features(:, indices);
    features = features';
    
    %Area under curve of SampEnMF (Attribute A according [5])
    results(:, 1) = trapz(features);

    %Area ratio of SampEnMF (Attribute Γ according [5])
    first = features;
    last = features;
    
    first = first( 1:floor(size(first,1)/2), :);
    last = last(((floor(size(last,1)/2)) + 1) : end, :);
    
    results(:, 2) = trapz(last) ./ trapz(first);
    
    %Skewness (Obliquity) of SampEnMF (Attribute O according [5])
    results(:, 3) = skewness(features);
 
    %Maximum point of SampEnMF (Attribute W according [5])
    [M,I]=max(features);
    results(:, 4) = M';
   
    %Maximum point scale of SampEnMF (Attribute E according [5])
    results(:, 5) = x(I);
    
end

