defmodule Wavelets.WaveletPacket do
  @moduledoc """
  Implementation of Wavelet Packet decomposition and reconstruction.
  """

  alias Wavelets.DWT
  alias Wavelets.Filter
  alias Wavelets.IDWT

  @doc """
  Performs 1D wavelet packet decomposition.
  Creates a tree structure where each node contains coefficients at different
  scales and frequency bands.
  """
  @spec decompose_1d(
          list(number),
          Filter.t(),
          non_neg_integer(),
          Wavelets.precision()
        ) ::
          map()
  def decompose_1d(signal, filter, levels, _precision \\ :double) do
    tree = %{
      {0, 0} => signal
    }

    do_decompose_1d(tree, filter, levels)
  end

  @doc """
  Reconstructs the signal from its wavelet packet decomposition.
  """
  @spec reconstruct_1d(map(), Filter.t(), Wavelets.precision()) :: list(number)
  def reconstruct_1d(tree, filter, _precision \\ :double) do
    max_level =
      tree
      |> Map.keys()
      |> Enum.map(fn {level, _} -> level end)
      |> Enum.max()

    do_reconstruct_1d(tree, filter, max_level)
  end

  # Private helper functions
  defp do_decompose_1d(tree, _filter, 0), do: tree

  defp do_decompose_1d(tree, filter, level) do
    current_level = level - 1
    max_nodes = trunc(:math.pow(2, current_level))

    new_nodes =
      for j <- 0..(max_nodes - 1),
          Map.has_key?(tree, {current_level, j}) do
        signal = Map.get(tree, {current_level, j})
        {approx, details} = DWT.forward_1d(signal, filter)

        %{
          {level, j * 2} => approx,
          {level, j * 2 + 1} => details
        }
      end

    if Enum.empty?(new_nodes) do
      tree
    else
      merged = Enum.reduce(new_nodes, %{}, &Map.merge(&2, &1))

      Map.merge(tree, merged)
      |> do_decompose_1d(filter, level - 1)
    end
  end

  defp do_reconstruct_1d(tree, _filter, 0) do
    Map.get(tree, {0, 0})
  end

  defp do_reconstruct_1d(tree, filter, level) when level > 0 do
    max_nodes = trunc(:math.pow(2, level - 1))

    new_nodes =
      for j <- 0..(max_nodes - 1),
          Map.has_key?(tree, {level, j * 2}) and
            Map.has_key?(tree, {level, j * 2 + 1}) do
        approx = Map.get(tree, {level, j * 2})
        details = Map.get(tree, {level, j * 2 + 1})

        # Ensure we have both coefficients before reconstruction
        if is_list(approx) and is_list(details) do
          reconstructed = IDWT.inverse_1d(approx, details, filter)
          {{level - 1, j}, reconstructed}
        end
      end
      |> Enum.reject(&is_nil/1)
      |> Enum.into(%{})

    tree = Map.merge(tree, new_nodes)
    do_reconstruct_1d(tree, filter, level - 1)
  end
end
