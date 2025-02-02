use std::{env, io};

use noodles_fasta::{self as fasta, record::Definition, repository::adapters::IndexedReader};

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let name = args.next().expect("missing name");

    let repository = IndexedReader::builder()
        .open(src)
        .map(fasta::Repository::new)?;

    let sequence = repository
        .get(&name)
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid name"))?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = fasta::Writer::new(handle);

    let record = fasta::Record::new(Definition::new(name, None), sequence);
    writer.write_record(&record)?;

    Ok(())
}
