use thiserror::Error;

#[derive(Error, Debug)]
pub enum MagicByteError {
    #[error("Insufficient bytes: got {got}, need at least {need}")]
    InsufficientBytes { got: usize, need: usize },
}

pub fn are_gzipped_magic_bytes(first_bytes: &[u8]) -> Result<bool, MagicByteError> {
    if first_bytes.len() < 2 {
        return Err(MagicByteError::InsufficientBytes {
            got: first_bytes.len(),
            need: 2,
        });
    }
    Ok(first_bytes[0] == 0x1f && first_bytes[1] == 0x8b)
}
