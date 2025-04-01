import tensorflow as tf
from keras import layers

class GraphDiscriminatorDiff(tf.keras.Model):
    def __init__(self, num_nodes, node_features):
        super().__init__()

        self.num_nodes = num_nodes
        self.node_features = node_features

        self.conv1 = layers.Dense(64, activation = 'relu')
        self.conv2 = layers.Dense(128, activation = 'relu')

        self.fc = layers.Dense(1, activation = 'sigmoid')

        self.optimizer = tf.keras.optimizers.SGD(learning_rate=0.0002, momentum=0.5)    
   
    def call(self, adj, node_features):
        x = self.conv1(node_features)
        x = self.conv2(x)
        output = self.fc(x)

        return output 
