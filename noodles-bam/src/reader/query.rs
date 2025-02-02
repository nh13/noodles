use std::io::{self, Read, Seek};

use noodles_bgzf::{self as bgzf, VirtualPosition};
use noodles_core::region::Interval;
use noodles_csi::index::reference_sequence::bin::Chunk;
use noodles_sam::alignment::Record;

use super::Reader;

enum State {
    Seek,
    Read(VirtualPosition),
    Done,
}

/// An iterator over records of a BAM reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'a, R>
where
    R: Read + Seek,
{
    reader: &'a mut Reader<bgzf::Reader<R>>,

    chunks: Vec<Chunk>,
    i: usize,

    reference_sequence_id: usize,
    interval: Interval,

    state: State,
    record: Record,
}

impl<'a, R> Query<'a, R>
where
    R: Read + Seek,
{
    pub(super) fn new(
        reader: &'a mut Reader<bgzf::Reader<R>>,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader,

            chunks,
            i: 0,

            reference_sequence_id,
            interval,

            state: State::Seek,
            record: Record::default(),
        }
    }

    fn next_record(&mut self) -> io::Result<Option<Record>> {
        self.reader.read_record(&mut self.record).map(|n| match n {
            0 => None,
            _ => Some(self.record.clone()),
        })
    }
}

impl<'a, R> Iterator for Query<'a, R>
where
    R: Read + Seek,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.state {
                State::Seek => {
                    self.state = match next_chunk(&self.chunks, &mut self.i) {
                        Some(chunk) => {
                            if let Err(e) = self.reader.seek(chunk.start()) {
                                return Some(Err(e));
                            }

                            State::Read(chunk.end())
                        }
                        None => State::Done,
                    }
                }
                State::Read(chunk_end) => match self.next_record() {
                    Ok(Some(record)) => {
                        if self.reader.virtual_position() >= chunk_end {
                            self.state = State::Seek;
                        }

                        if intersects(&record, self.reference_sequence_id, self.interval) {
                            return Some(Ok(record));
                        }
                    }
                    Ok(None) => self.state = State::Seek,
                    Err(e) => return Some(Err(e)),
                },
                State::Done => return None,
            }
        }
    }
}

pub(crate) fn next_chunk(chunks: &[Chunk], i: &mut usize) -> Option<Chunk> {
    let chunk = chunks.get(*i).copied();
    *i += 1;
    chunk
}

pub(crate) fn intersects(
    record: &Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> bool {
    match (
        record.reference_sequence_id(),
        record.alignment_start(),
        record.alignment_end(),
    ) {
        (Some(id), Some(start), Some(end)) => {
            let alignment_interval = (start..=end).into();
            id == reference_sequence_id && region_interval.intersects(alignment_interval)
        }
        _ => false,
    }
}
