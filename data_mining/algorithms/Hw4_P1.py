
# coding: utf-8

# In[170]:


## SVM via Gradient Descent 
## Implemented in Jupyter notebooks

import numpy as np
import matplotlib.pyplot as plt
import timeit


# In[174]:


def Batch_GD(data, target):

    delta_cost = 1
    eta = 0.0000003
    epi = 0.25
    k = 0

    w = np.zeros(len(data[0,:]))
    w_new = w
    b = 0
    C = 100

    f_k = 0.5*np.sum(w**2) + C*np.sum(np.maximum(0,1-target*(np.dot(data,w) + b)))
    cost_func = [f_k]
    
    while delta_cost > epi:
        
        for j in range(len(w)):
            w_new[j] = w[j] - eta * (w[j] + C*np.sum(((target*(np.dot(data,w)+b))<1)*(-target*data[:,j])))
            
        b = b - eta * C*np.sum(((target*(np.dot(data,w)+b))<1)*-target)
        k = k+1
        
        f_k_1 = f_k
        f_k = 0.5*np.sum(w_new**2) + C*np.sum(np.maximum(0,1-target*(np.dot(data,w_new) + b)))
        
        delta_cost = abs(f_k_1-f_k)*100/f_k_1
        
        w = w_new
        cost_func.append(f_k)
        
    return k, cost_func      


# In[137]:


def Stochastic_GD(data, target):

    delta_cost = 0
    eta = 0.0001
    epi = 0.001
    i = 1
    k = 0
    n = len(data[:,0])

    w = np.zeros(len(data[0,:]))
    w_new = w
    b = 0
    C = 100
    
    indices = np.arange(target.shape[0])
    np.random.shuffle(indices)
    data = data[indices]
    target = target[indices]
    
    f_k = 0.5*np.sum(w**2) + C*np.sum(np.maximum(0,1-target*(np.dot(data,w) + b)))
    cost_func = [f_k]
    
    while delta_cost > epi or k < 1:
        
        for j in range(len(w)):
            w_new[j] = w[j] - eta * (w[j] + C*((target[i]*(np.dot(data[i,:],w)+b))<1)*(-target[i]*data[i,j]))
            
        b = b - eta * C*(((target[i]*(np.dot(data[i,:],w)+b))<1)*-target[i])
        i = i % n + 1
        k = k+1
        
        f_k_1 = f_k
        f_k = 0.5*np.sum(w_new**2) + C*np.sum(np.maximum(0,1-target*(np.dot(data,w_new) + b)))
        
        delta_cost = 0.5 * delta_cost + 0.5 * abs(f_k_1-f_k)*100/f_k_1
        
        w = w_new
        cost_func.append(f_k)
        
    return k, cost_func   


# In[158]:


def Mini_Batch_GD(data, target):

    delta_cost = 0
    eta = 0.00001
    epi = 0.01
    batch_size = 20
    k = 0
    n = len(data[:,0])
    l = 0

    w = np.zeros(len(data[0,:]))
    w_new = w
    b = 0
    C = 100
    
    indices = np.arange(target.shape[0])
    np.random.shuffle(indices)
    data = data[indices]
    target = target[indices]
    
    f_k = 0.5*np.sum(w**2) + C*np.sum(np.maximum(0,1-target*(np.dot(data,w) + b)))
    cost_func = [f_k]
    
    while delta_cost > epi or k < 1:
        
        i = int(l*batch_size+1)
        end = int(min(n,(l+1)*batch_size))
        
        for j in range(len(w)):
            w_new[j] = w[j] - eta * (w[j] + C*np.sum(((target[i:end]*(np.dot(data[i:end,:],w)+b))<1)
                                            *(-target[i:end]*data[i:end,j])))
            
        b = b - eta * C*np.sum(((target[i:end]*(np.dot(data[i:end,:],w)+b))<1)*-target[i:end])
        l = (l+1) % ((n+batch_size - 1)/batch_size)
        k = k+1
        
        f_k_1 = f_k
        f_k = 0.5*np.sum(w_new**2) + C*np.sum(np.maximum(0,1-target*(np.dot(data,w_new) + b)))
        
        delta_cost = 0.5 * delta_cost + 0.5 * abs(f_k_1-f_k)*100/f_k_1
        
        w = w_new
        cost_func.append(f_k)
        
    return k, cost_func  


# In[173]:


data = np.loadtxt("features.txt", delimiter= ",")
target = np.loadtxt("target.txt")


# In[175]:


start = timeit.default_timer()
k, cost_func = Batch_GD(data, target)
stop = timeit.default_timer()

BGD_runtime = stop-start


# In[176]:


BGD_runtime


# In[194]:


start = timeit.default_timer()
SGD_k, SGD_cost_func = Stochastic_GD(data, target)
stop = timeit.default_timer()

SGD_runtime = stop-start
plt.plot(SGD_cost_func)


# In[195]:


SGD_runtime



# In[179]:


start = timeit.default_timer()
MBGD_k, MBGD_cost_func = Mini_Batch_GD(data, target)
stop = timeit.default_timer()

MBGD_runtime = stop-start
plt.plot(MBGD_cost_func)


# In[180]:


MBGD_runtime


# In[168]:


plt.figure()
plt.plot(cost_func, label="Batch GD")
plt.plot(SGD_cost_func, label="Stochastic GD")
plt.plot(MBGD_cost_func, label="Mini Batch GD")
plt.legend()
plt.xlabel("Number of Iterations (k)")
plt.ylabel("Cost Function")
plt.title("Implmentation of SVM with \n Various Gradient Descent Methods")

