pub(crate) mod compression_header;
mod slice;

use std::{
    cmp,
    io::{self, Write},
};

use self::compression_header::write_compression_header;
use crate::{
    container::Block,
    data_container::{Header, ReferenceSequenceContext, Slice},
    DataContainer,
};

pub fn write_data_container<W>(
    writer: &mut W,
    data_container: &DataContainer,
    base_count: u64,
) -> io::Result<()>
where
    W: Write,
{
    use super::container::{write_block, write_header};

    let (header, blocks) = build_container(data_container, base_count)?;

    write_header(writer, &header)?;

    for block in blocks {
        write_block(writer, &block)?;
    }

    Ok(())
}

fn build_container(
    data_container: &DataContainer,
    base_count: u64,
) -> io::Result<(Header, Vec<Block>)> {
    use crate::container::block::ContentType;

    let mut buf = Vec::new();
    write_compression_header(&mut buf, data_container.compression_header())?;

    let block = Block::builder()
        .set_content_type(ContentType::CompressionHeader)
        .set_uncompressed_len(buf.len())
        .set_data(buf.into())
        .build();

    let mut blocks = vec![block];
    let mut landmarks = Vec::new();

    let container_reference_sequence_context =
        build_container_reference_sequence_context(data_container.slices())?;

    let mut container_record_count = 0;
    let container_record_counter = data_container
        .slices()
        .first()
        .map(|s| s.header().record_counter())
        .expect("no slices in builder");

    for slice in data_container.slices() {
        let slice_header = slice.header();

        container_record_count += slice_header.record_count() as i32;

        let mut slice_len = 0;

        let mut slice_header_buf = Vec::new();
        self::slice::write_header(&mut slice_header_buf, slice.header())?;

        let slice_header_block = Block::builder()
            .set_content_type(ContentType::SliceHeader)
            .set_uncompressed_len(slice_header_buf.len())
            .set_data(slice_header_buf.into())
            .build();

        slice_len += slice_header_block.len();
        blocks.push(slice_header_block);

        blocks.push(slice.core_data_block().clone());
        slice_len += slice.core_data_block().len();

        for external_block in slice.external_blocks() {
            blocks.push(external_block.clone());
            slice_len += external_block.len();
        }

        let last_landmark = landmarks.last().copied().unwrap_or(0);
        let landmark = last_landmark + slice_len;
        landmarks.push(landmark);
    }

    let len = blocks.iter().map(|b| b.len()).sum();

    let header = Header::builder()
        .set_length(len)
        .set_reference_sequence_context(container_reference_sequence_context)
        .set_record_count(container_record_count)
        .set_record_counter(container_record_counter)
        .set_base_count(base_count)
        .set_block_count(blocks.len())
        .set_landmarks(landmarks)
        .build();

    Ok((header, blocks))
}

fn build_container_reference_sequence_context(
    slices: &[Slice],
) -> io::Result<ReferenceSequenceContext> {
    assert!(!slices.is_empty());

    let first_slice = slices.first().expect("slices cannot be empty");
    let mut container_reference_sequence_context =
        first_slice.header().reference_sequence_context();

    for slice in slices.iter().skip(1) {
        let slice_reference_sequence_context = slice.header().reference_sequence_context();

        match (
            container_reference_sequence_context,
            slice_reference_sequence_context,
        ) {
            (
                ReferenceSequenceContext::Some(container_context),
                ReferenceSequenceContext::Some(slice_context),
            ) if container_context.reference_sequence_id()
                == slice_context.reference_sequence_id() =>
            {
                let alignment_start = cmp::min(
                    container_context.alignment_start(),
                    slice_context.alignment_start(),
                );

                let alignment_end = cmp::max(
                    container_context.alignment_end(),
                    slice_context.alignment_end(),
                );

                container_reference_sequence_context = ReferenceSequenceContext::some(
                    container_context.reference_sequence_id(),
                    alignment_start,
                    alignment_end,
                );
            }
            (ReferenceSequenceContext::None, ReferenceSequenceContext::None) => {}
            (ReferenceSequenceContext::Many, ReferenceSequenceContext::Many) => {}
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "invalid slice reference sequence context: expected {:?}, got {:?}",
                        container_reference_sequence_context, slice_reference_sequence_context
                    ),
                ));
            }
        }
    }

    Ok(container_reference_sequence_context)
}
