function [abundance_matrix, mixed_data_with_grain] = simulate_Mixed_abundance_data(normalizedMatrix, endmember_matrix, grain_size, min_abundance, max_abundance)
   % Scale to desired range [a,b]
   % First, determine current range after normalization
   currentMin = min(normalizedMatrix, [], 'all');
   currentMax = max(normalizedMatrix, [], 'all');

   % Apply linear scaling
   abundance_matrix = min_abundance + (max_abundance - min_abundance) * (normalizedMatrix - currentMin) / (currentMax - currentMin);

   % Renormalize to ensure rows still sum to 1
   rowSumsFinal = sum(abundance_matrix, 2);
   abundance_matrix = abundance_matrix ./ rowSumsFinal;
   %% Generate Mixed Distributions for Each Specimen
   mixed_distributions = abundance_matrix * endmember_matrix;
   % Method 1: Save as CSV with grain size as first column
   mixed_data_with_grain = [grain_size, mixed_distributions'];
end

