%DEMO Demonstration of parametric t-SNE


    % Load MNIST dataset
    load 'mnist_train.mat'
    load 'mnist_test.mat'
    
    % Set perplexity and network structure
    perplexity = 30;
    %layers = [500 500 2000 2];
    layers = [50 50 200 2];
    numSamples = 20000;
    % Train the parametric t-SNE network
    %[network, err] = train_par_tsne(train_X, train_labels, test_X, test_labels, layers, 'CD1', perplexity);
    [network] = train_par_tsne(train_X(1:numSamples,:), layers, 'CD1', perplexity);
    
    % Construct training and test embeddings
    mapped_train_X = run_data_through_network(network, train_X);
    mapped_test_X  = run_data_through_network(network, test_X);
    
    % Compute 1-NN error and trustworthiness
    %disp(['1-NN error: ' num2str(knn_error(mapped_train_X, train_labels, mapped_test_X, test_labels, 1))]);
    %disp(['Trustworthiness: ' num2str(trustworthiness(test_X, mapped_test_X, 12))]);
    
    %% Plot test embedding
    figure(1);
    subplot(2,1,1);
    scatter(mapped_train_X(1:numSamples,1), mapped_train_X(1:numSamples,2), 9, train_labels(1:numSamples));    
    subplot(2,1,2);
    scatter(mapped_test_X(1:numSamples,1), mapped_test_X(1:numSamples,2), 9, test_labels(1:numSamples));
    title('Embedding of test data');