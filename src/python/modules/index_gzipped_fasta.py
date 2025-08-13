#!/usr/bin/env python3

"""
Creates a binary index for gzip-compressed exon Fasta file
"""

from typing import List, TextIO
from shared import CONTEXT_SETTINGS

import click
import gzip
import struct
import os
# import zlib

HEADER_START: bytes = b'>'
NEWLINE: bytes = b'\n'
QUERY_EXON: str = b'query_exon'

WRITE_COUNTER: int = 0

def write_index_line(
    proj: bytes, exon_num: int, start_byte: int, end_byte: int, file: TextIO
) -> None:
    global WRITE_COUNTER
    proj_byte_len: int = len(proj)
    file.write(struct.pack('>I', proj_byte_len))
    file.write(proj)
    file.write(struct.pack('>Q', exon_num))
    file.write(struct.pack('>Q', start_byte))
    file.write(struct.pack('>Q', end_byte))
    WRITE_COUNTER += 1


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'input_',
    type=click.Path(exists=True),
    metavar='INPUT'
)

def main(input_: click.Path) -> None:
    global WRITE_COUNTER
    in_dir: str = os.path.dirname(input_)
    out_file: str = '.' + os.path.basename(input_) + '.ix'
    output_: str = os.path.join(in_dir, out_file)
    with gzip.open(input_, 'rb') as ih, open(output_, 'wb') as oh:
        # decompressor = zlib.decompressobj(32 + zlib.MAX_WBITS)
        buffer: bytes = b''
        proj: bytes = b''
        exon_num: int = ''
        start_byte: int = 0
        current_position: int = 0
        processed_len: int = 0
        while True:
            # if WRITE_COUNTER and not WRITE_COUNTER % 1000:
            #     print(f'{WRITE_COUNTER} records already indexed')
            chunk: bytes = ih.read(1024)
            if not chunk:
                break
            buffer += chunk#decompressed_data
            while True:
                newline_index: int = buffer.find(NEWLINE)
                if newline_index < 0:
                    break
                # print(f'{buffer=}')
                line: bytes = buffer[:newline_index + 1]
                buffer = buffer[newline_index + 1:]
                # print(f'{line=}')
                # print(f'{buffer=}')
                # current_position = ih.tell() - len(buffer) + len(line) + processed_len
                current_position: int = ih.tell() - len(buffer) - newline_index - 1 # + processed_len
                # print(f'{start_byte=}, {current_position=}, {ih.tell()=}, {len(buffer)=}, {processed_len=}')
                if line.startswith(b'>'):
                    if proj:
                        write_index_line(proj, exon_num, start_byte, current_position, oh)
                        proj = b''
                        exon_num = 0
                    header_split: List[bytes] = line.rstrip().split(b' | ')
                    # print(f'{header_split=}')
                    source: bytes = header_split[-1]
                    if source == QUERY_EXON:
                        # print('Updating starting byte')
                        proj: bytes = header_split[0][1:]
                        exon_num = int(header_split[1])
                        start_byte = current_position
            processed_len += len(chunk) - len(buffer)
            # print(f'Finishing chunk (residual buffer length = {len(buffer)})\n')
        #     chunk: bytes = ih.read(1024)
        #     if not chunk:
        #         break
        #     decompressed_data: bytes = decompressor.decompress(chunk)
        #     buffer += decompressed_data
        #     while True:
        #         newline_index: int = buffer.find(NEWLINE)
        #         if newline_index < 0:
        #             break
        #         line: bytes = buffer[:newline_index + 1]
        #         buffer = buffer[newline_index + 1:]
        #         current_position = ih.tell() - len(chunk) + len(decompressed_data) - len(buffer)
        #         if line.startswith(b'>'):
        #             if proj:
        #                 write_index_line(proj, exon_num, start_byte, current_position, oh)
        #                 proj = ''
        #                 exon_num = 0
        #             header_split: List[bytes] = line.rstrip().split(b' | ')
        #             source: bytes = header_split[-1]
        #             if source == QUERY_EXON:
        #                 proj: bytes = header_split[0][1:]
        #                 exon_num = int(header_split[1])
        #                 start_byte = current_position
        # with gzip.GzipFile(fileobj=ih) as gz:
        #     proj: str = ''
        #     exon_num: int = ''
        #     start_byte: int = 0
        #     while True:
        #         if not WRITE_COUNTER % 1000:
        #             print(f'{WRITE_COUNTER} records already indexed')
        #         line_start_pos: int = ih.tell()
        #         print(f'{line_start_pos=}')
        #         line: str = gz.readline()
        #         if not line:
        #             if proj:
        #                 write_index_line(proj, exon_num, start_byte, line_start_pos, oh)
        #             break
        #         if line.startswith(b'>'):
        #             if proj:
        #                 write_index_line(proj, exon_num, start_byte, line_start_pos, oh)
        #                 proj = ''
        #                 exon_num = 0
        #             header_split: List[bytes] = line.rstrip().split(b' | ')
        #             source: bytes = header_split[-1]
        #             if source == QUERY_EXON:
        #                 proj: bytes = header_split[0][1:]
        #                 exon_num = int(header_split[1])
        #                 start_byte = line_start_pos


if __name__ == '__main__':
    main()