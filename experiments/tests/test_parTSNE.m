% Set perplexity and network structure
train_X = output{1}.optima;
test_X = output{1}.conceptSelection.members;

% Train the parametric t-SNE network
perplexity = 50;
layers = [10 10 50 2];
[network] = train_par_tsne(output{1}.optima, layers, 'CD1', perplexity);

% Construct training and test embeddings
mapped_train_X = run_data_through_network(network, train_X);
mapped_test_X  = run_data_through_network(network, test_X);

disp(['Trustworthiness: ' num2str(trustworthiness(test_X, mapped_test_X, 12))]);

%% Plot test embedding
figure(1);
hold off;
scatter(mapped_train_X(:,1), mapped_train_X(:,2), 18, 'k', 'filled');
hold on;
scatter(mapped_test_X(:,1), mapped_test_X(:,2), 9, 'r', 'filled');
title('Embedding of test data');