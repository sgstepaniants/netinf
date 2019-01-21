import shelve

#******************************************************************************
# Read Neural Net Results
#******************************************************************************
filepath = 'data/free/n=2/p=0.75/'
results = shelve.open(filepath + 'results')

exclude = ['learning_rate', 'weight_decay', 'train_err', 'test_err']
for key in sorted(results):
    if key not in exclude:
        print('%s: %s' % (key, results[key]))
results.close()
