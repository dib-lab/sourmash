use std::io::Write;

use crate::Error;

pub trait Update<O> {
    fn update(&self, other: &mut O) -> Result<(), Error>;
}

pub trait ToWriter {
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write;
}
