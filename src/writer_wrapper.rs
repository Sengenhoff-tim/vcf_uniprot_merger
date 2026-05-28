use flate2::write::GzEncoder;
use std::io::Write;

pub enum WriterWrapper<W: Write> {
    Compressed(GzEncoder<W>),
    Uncompressed(W),
}

impl<W: Write> Write for WriterWrapper<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            WriterWrapper::Compressed(enc) => enc.write(buf),
            WriterWrapper::Uncompressed(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            WriterWrapper::Compressed(enc) => enc.flush(),
            WriterWrapper::Uncompressed(w) => w.flush(),
        }
    }
}
