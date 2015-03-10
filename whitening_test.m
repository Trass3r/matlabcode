X = rand(100,20); % 100 instance with 20 features

figure;
imagesc(cov(X)); colorbar; title('original covariance matrix');

[Xwh, mu, invMat, whMat] = whiten(X,0.0001);

figure;
imagesc(cov(Xwh)); colorbar; title('whitened covariance matrix');

Xwh2 = (X-repmat(mean(X), size(X,1),1))*whMat; 
figure;
plot(sum((Xwh-Xwh2).^2),'-rx'); title('reconstructed whitening error (should be 0)');

Xrec = Xwh*invMat + repmat(mu, size(X,1),1);
figure;
plot(sum((X-Xrec).^2),'-rx'); ylim([-1,1]); title('reconstructed data error (should be zero)');
