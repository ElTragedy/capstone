import os

def save_checkpoint(generator, discriminator, epoch, output_dir):
    gen_path = os.path.join(output_dir, f'generator_epoch_{epoch}')
    disc_path = os.path.join(output_dir, f'discriminator_epoch_{epoch}')
    
    os.makedirs(gen_path, exist_ok=True)
    os.makedirs(disc_path, exist_ok=True)
    
    generator.save_weights(os.path.join(gen_path, 'weights.weights.h5'))
    discriminator.save_weights(os.path.join(disc_path, 'weights.weights.h5'))
    
    print(f"Checkpoint saved for epoch {epoch} to {output_dir}")
