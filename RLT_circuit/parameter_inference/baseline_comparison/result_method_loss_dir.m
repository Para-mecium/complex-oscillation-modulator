function resultDir = result_method_loss_dir(resultsRoot, methodName, lossName)
methodToken = result_name_token(methodName);
lossToken = result_name_token(lossName);
resultDir = fullfile(resultsRoot, methodToken, lossToken);
end
