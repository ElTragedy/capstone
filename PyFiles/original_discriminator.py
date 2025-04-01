import tensorflow as tf
<<<<<<< HEAD
from keras import layers

#base graph discriminator without pooling
class GraphDiscriminatorNon(tf.keras.Model):
    def __init__(self, num_nodes, node_features):
        super().__init__()
=======
from tensorflow.keras import layers, models

class GraphDiscriminator(tf.keras.Model):
    def __init__(self, num_nodes, node_features):
        super().__init__()
        self.input_dim = num_nodes * num_nodes
        self.condition_dim = num_nodes * node_features
>>>>>>> 3d65a8acb99b7207715346f054cda2f2e3f95ce8

        self.mlp = models.Sequential([
            layers.Dense(256, input_shape = (self.input_dim + self.condition_dim,), activation='leaky_relu'),
            layers.Dense(128, activation='leaky_relu'),
            #layers.BatchNormalization(),
            layers.Dense(1, activation='sigmoid')
        ])

        self.optimizer = tf.keras.optimizers.SGD(learning_rate=0.000099, momentum=0.5)

<<<<<<< HEAD
        #fully connected output layer (sigmoid for binary classification)
        self.fc = layers.Dense(1, activation = 'sigmoid')

        #optmizer for training the discriminator
        self.optimizer = tf.keras.optimizers.SGD(learning_rate=0.0002, beta_1=0.5)
        
    #forward pass for discriminator; adg - adjacency matrix, node_features: node feature matrix
    def call(self, adj, node_features):
        x = self.conv1(node_features)
        x = self.conv2(x)
        output = self.fc(x)
        return output #probability score indicating real or fake molecular graph
=======

    def call(self, adj, node_features):
        batch_size = tf.shape(adj)[0]

        adj_flat = tf.reshape(adj, (batch_size, self.input_dim))
        node_features_flat = tf.reshape(node_features, (batch_size, self.condition_dim))

        combined = tf.concat([adj_flat, node_features_flat], axis = 1)

        validity = self.mlp(combined)

        return validity
>>>>>>> 3d65a8acb99b7207715346f054cda2f2e3f95ce8
