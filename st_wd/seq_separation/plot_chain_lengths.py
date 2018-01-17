import matplotlib.pyplot as plt
import pickle as pkl
plt.style.use('ggplot')

ls = pkl.load(open('db_chain_lengths.pkl','rb'))
plt.hist(ls, bins=100)
plt.xlim(0, 1200)
plt.xticks([i*20 for i in range(60)])
plt.xticks(rotation=90)
plt.title('chain lengths in nrPDB')
plt.xlabel('chain length')
plt.ylabel('counts')
plt.tight_layout()
plt.show()
