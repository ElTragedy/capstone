import tensorflow as tf
from keras import layers, models

#base graph discriminator without pooling
class GraphDiscriminator(tf.keras.Model):
    def __init__(self, num_nodes, node_features):
        super(GraphDiscriminator, self).__init__()

        self.num_nodes = num_nodes
        self.node_features = node_features

        #dense layers to process node features
        self.conv1 = layers.Dense(64, activation = 'relu')
        self.conv2 = layers.Dense(128, activation = 'relu')

        #global pooling to combine features across nodes
        self.pool = layers.GlobalAveragePooling1D()

        #fully connected output layer (sigmoid for binary classification)
        self.fc = layers.Dense(1, activation = 'sigmoid')

        #optmizer for training the discriminator
        self.optimizer = tf.keras.optimizers.Adam(learning_rate=0.0002, beta_1=0.5)
        
    #forward pass for discriminator; adg - adjacency matrix, node_features: node feature matrix
    def call(self, adj, node_features):
        x = self.conv1(node_features)
        x = self.conv2(x)
        x = self.pool(x)
        output = self.fc(x)

        return output #probability score indicating real or fake molecular graph
