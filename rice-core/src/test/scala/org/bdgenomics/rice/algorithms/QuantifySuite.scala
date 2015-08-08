/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.bdgenomics.rice.algorithms

import org.apache.spark.rdd.RDD
import org.bdgenomics.formats.avro.{ Contig, NucleotideContigFragment, AlignmentRecord }
import org.bdgenomics.adam.models.{ Exon, ReferenceRegion, Transcript, CDS, UTR }
import org.bdgenomics.rice.models.{ KmerIndex, IndexMap }
import org.bdgenomics.rice.algorithms.alignment.AlignmentModel
import org.bdgenomics.rice.utils.riceFunSuite
import net.fnothaft.ananas.models._
import net.fnothaft.ananas.avro.{Kmer, Backing}
import net.fnothaft.ananas.debruijn.ColoredDeBruijnGraph

object TestAlignmentModel extends AlignmentModel {

  def processRead(iter: Iterator[CanonicalKmer],
                  kmerIndex: KmerIndex): Map[String, Double] = {
      
    iter.toList
        .flatMap( c => kmerIndex.getTranscripts(c) )
        .groupBy(_._1) // Map[ TranscriptID -> List[(TranscriptId, Occurrences)] ] 
        .map( v => v._2.reduce( (a, b) => (a._1, a._2 + b._2) ) )
        .map( v => (v._1, v._2.toDouble) )
  }  
}

class QuantifySuite extends riceFunSuite {

  // Blatant copy of function from ananas/models/ContigFragment to get around protected status
  def buildFromNCF(fragment: NucleotideContigFragment): ContigFragment = {
    // is this the last fragment in a contig?
    val isLast = fragment.getFragmentNumber == (fragment.getNumberOfFragmentsInContig - 1)

    val sequence = IntMer.fromSequence(fragment.getFragmentSequence)
      .map(_.asInstanceOf[CanonicalKmer])

    new ContigFragment(fragment.getContig.getContigName,
                       sequence,
                       isLast,
                       Option(fragment.getFragmentStartPosition).fold(0)(_.toInt))
  }

  def createTestIndex(sequence: String = "ACACTGTGGGTACACTACGAGA") : (Map[Long, Map[String, Long]], Map[String, Transcript]) = {
    val ncf = NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
      .setContigName("ctg")
      .build())
      .setFragmentNumber(0)
      .setNumberOfFragmentsInContig(1)
      .setFragmentSequence(sequence)
      .build()

    val frag = sc.parallelize( Seq( buildFromNCF(ncf) ) ) 

    val tx = Seq( Transcript("one", Seq("one"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]()) ,
                  Transcript("two", Seq("two"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]()) )
    val transcripts = sc.parallelize(tx)

    Index(frag, transcripts)
  }

  sparkTest("Simple Test of Index") {
    val testSeq = "ACACTGTGGGTACACTACGAGA"
    val (imap, tmap) = createTestIndex(testSeq)

    // Test kmer mapping
    val imers = testSeq.sliding(16).map(s => IntMer(s))
 
    assert( imap.size == 7 ) // 7 kmers of length 16
    assert( imers.forall(i => imap(i.longHash)("ctg") == 1) )

    // Test transcript mapping
    assert( tmap("one").id == "one")
    assert( tmap("two").id == "two")
  }

  sparkTest("Simple Test of Mapper") {
    val testSeq = "ACACTGTGGGTACACTACGAGA"
    val ar = Array({AlignmentRecord.newBuilder()
                            .setSequence(testSeq)
                            .build()})

    val reads = sc.parallelize(ar)

    val (imap, tmap) = createTestIndex(testSeq)

    val kmerIndex = IndexMap(16, imap)

    val m = Mapper(reads, kmerIndex, TestAlignmentModel).collect()

    // Only one read, so only one item in the map
    assert( m.size == 1 )

    // The one item should map to ( "ctg" -> number of occurrences of all kmers in contig ctg = 7, as all 7 kmers are from the read)
    assert( m.forall( v => v._2("ctg") == 7.toDouble ) )

  }

  sparkTest("Less Simple Test of Mapper") {
    val testSeq1 = "ACACTGTGGGTACACTACGAGA"
    val testSeq2 = "CCAGTGACTGGAAAAA"
    val ar = Array({AlignmentRecord.newBuilder()
                            .setSequence(testSeq1)
                            .build()},
                    {AlignmentRecord.newBuilder()
                            .setSequence(testSeq2)
                            .build()})

    val reads = sc.parallelize(ar)

    val (imap, tmap) = createTestIndex(testSeq1 + testSeq2)

    val kmerIndex = IndexMap(16, imap)

    val m = Mapper(reads, kmerIndex, TestAlignmentModel).collect()

    // Two reads, so two items in the map
    assert( m.size == 2 )

    // Assert that readIDs are unique
    assert( m(0)._1 != m(1)._1 )

    // Should have either 7 or 1 kmers per read (all kmers are from the same transcript "ctg")
    assert( m.forall( v => {v._2("ctg") == 7.toDouble || v._2("ctg") == 1} ) )

    // Should record 8 kmers in total
    assert( m(0)._2("ctg") + m(1)._2("ctg") == 8)

  }
}
